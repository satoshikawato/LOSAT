//! The even-gap linking of the HSPs of one subject with sum statistics (NCBI core/link_hsps.c
//! `BLAST_LinkHsps` with `longest_intron <= 0`, `s_BlastEvenGapLinkHSPs`) for an untranslated
//! subject: psiblast's `-sum_stats true` (a gapped search, untranslated query) and BLASTX's
//! `-ungapped` search with a `-max_intron_length` of at most 4 (an ungapped search, translated
//! query: `longest_intron` is then 0, blast_parameters.c:809-814).
//!
//! NCBI's blastp and blastn have no `-sum_stats`, and the other linking searches take other forms:
//! TBLASTX links ungapped HSPs of translated sequences (`algorithm::tblastx::sum_stats_linking`, its
//! own kernels for its `UngappedHit`), and the gapped BLASTX and TBLASTN searches use the
//! uneven-gap linking (`s_BlastUnevenGapLinkHSPs`). The cutoffs of the ungapped searches' linking
//! (`CalculateLinkHSPCutoffs`) are `calculate_link_hsp_cutoffs`'s, which TBLASTX calls too.
//!
//! The port is literal, with NCBI's shortcuts (`use_current_max`, `path_changed`, the `changed`
//! flags, `next_larger`): they decide which of equal sums a chain takes, so the result depends on
//! them. The sums are `i64` where NCBI keeps `Int4` (a sum past `Int4` is undefined behaviour in
//! NCBI; Owner decision 2026-10-10). The callers give the HSPs in the list's order and read back the
//! list NCBI builds (`LinkedHsp`, in NCBI's order of `hsp_list->hsp_array` after the linking), then
//! sort it by score (`Blast_HSPListSortByScore`) as `BLAST_LinkHsps` does.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1760-1814

use crate::core::blast_util::ncbi_int4_from_double;
use crate::utils::glibc_fma_libm::GlibcLibm;
use std::cmp::Ordering;

use crate::stats::sum_statistics::defaults::{GAP_SIZE, OVERLAP_SIZE};
use crate::stats::sum_statistics::{gap_decay_divisor, ncbi_large_gap_sum_e, small_gap_sum_e};
use crate::stats::KarlinParams;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:66-71
// With a gap probability of 1, the preliminary stage (which also links with small gaps) never
// takes a large-gap set of more than one HSP (`1 - gap_prob == 0` gives it `INT4_MAX`); the
// traceback links with large gaps only and does not divide by the probability. `BLAST_GAP_SIZE`
// and `BLAST_OVERLAP_SIZE` are `stats::sum_statistics::defaults`'.
/// `BLAST_GAP_PROB_GAPPED`.
pub(crate) const BLAST_GAP_PROB_GAPPED: f64 = 1.0;
/// `BLAST_GAP_PROB`, the gap probability of an ungapped search.
const BLAST_GAP_PROB: f64 = 0.5;
/// `BLAST_GAP_DECAY_RATE_GAPPED`.
pub(crate) const BLAST_GAP_DECAY_RATE_GAPPED: f64 = 0.1;
/// `BLAST_GAP_DECAY_RATE`, the gap decay rate of an ungapped search.
pub(crate) const BLAST_GAP_DECAY_RATE: f64 = 0.5;

/// The largest gap of the uneven-gap linking (`BlastLinkHSPParameters::longest_intron`, in letters
/// of the translated sequence) that a search with sum statistics gets from `-max_intron_length`
/// (`BlastHitSavingOptions::longest_intron`, in nucleotides), for the translated searches other
/// than TBLASTX (BLASTX, TBLASTN and its `-in_pssm` search). `None`: NCBI frees the linking
/// parameters and clears `params->do_sum_stats`, so the search runs as it does without sum
/// statistics (a gapped search with a value of 1 to 4). 0 (an ungapped search with a value up to
/// 4): the even-gap linking (`BLAST_LinkHsps`, link_hsps.c:1779 `longest_intron <= 0`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:77-78
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:788-815
/// The option is not negative (the command line refuses a negative value, as NCBI's does), so
/// `longest_intron - 2` cannot pass the range of `Int4`; C's division truncates toward zero, as
/// Rust's does.
pub(crate) fn translated_longest_intron(max_intron_length: i32, gapped: bool) -> Option<i32> {
    const DEFAULT_LONGEST_INTRON: i32 = 122;
    let max_protein_gap = (max_intron_length - 2) / 3;
    if !gapped {
        Some(max_protein_gap.max(0))
    } else if max_intron_length == 0 {
        Some((DEFAULT_LONGEST_INTRON - 2) / 3)
    } else if max_protein_gap <= 0 {
        None
    } else {
        Some(max_protein_gap)
    }
}

/// The linking parameters (`BlastLinkHSPParameters`) of a protein search with sum statistics and
/// even-gap linking: a gapped one (`preliminary`, `traceback`; psiblast) or an ungapped one
/// (`ungapped`; BLASTX), and whether the query is translated.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:609-618
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:633-637
/// A gapped search sets no `cutoff_big_gap` (`CalculateLinkHSPCutoffs` runs only for ungapped
/// searches, blast_engine.c:1448-1455), so it stays at the 0 of the `calloc`. The preliminary stage
/// sets `cutoff_small_gap` to the smallest cutoff of the initial words
/// (`BlastLinkHSPParametersUpdate`, blast_engine.c:1399, blast_setup.c:1024); the traceback and
/// the composition-based redo link with the parameters of their own setup
/// (`SThreadLocalDataArraySetup` -> `BLAST_GapAlignSetUp` -> `BlastHitSavingParametersNew`), which
/// no update reaches, so their `cutoff_small_gap` is 0 too (`EvenGapLinkParams::traceback`).
/// `gap_size` and `overlap_size` are `BLAST_GAP_SIZE` and `BLAST_OVERLAP_SIZE` in every search (no
/// option or later setup changes them), which `even_gap_link_hsps` reads as the constants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EvenGapLinkParams {
    /// `gap_prob`.
    pub(crate) gap_prob: f64,
    /// `gap_decay_rate`.
    pub(crate) gap_decay_rate: f64,
    /// `cutoff_small_gap`: 0 links with large gaps only (`ignore_small_gaps`).
    pub(crate) cutoff_small_gap: i32,
    /// `cutoff_big_gap`.
    pub(crate) cutoff_big_gap: i32,
    /// `kTranslatedQuery` (`Blast_QueryIsTranslated`): the HSPs of the three frames of a strand
    /// (contexts `3k..3k+3`) form one query frame of the linking (`strand_factor` 3), sorted by
    /// `s_RevCompareHSPsTbx` and at the end by `s_FwdCompareHSPsTransl` (`even_gap_link_hsps`
    /// says why NCBI's `s_RevCompareHSPsTransl` before it changes nothing).
    pub(crate) query_translated: bool,
    /// `Blast_SubjectIsTranslated` (TBLASTN's ungapped search): the subject length passed to the
    /// linking is in nucleotides, and each query frame divides it and the context's length
    /// adjustment by `CODON_LENGTH` (`link_frame`).
    pub(crate) subject_translated: bool,
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
            cutoff_small_gap: 0,
            cutoff_big_gap: 0,
            query_translated: false,
            subject_translated: false,
        }
    }

    /// The parameters of an ungapped search's linking with one subject: the cutoffs and the gap
    /// probability `CalculateLinkHSPCutoffs` sets (`calculate_link_hsp_cutoffs`) on the parameters
    /// of `BlastLinkHSPParametersNew` for an ungapped search (`gap_decay_rate`
    /// `BLAST_GAP_DECAY_RATE`, blast_parameters.c:609-617).
    pub(crate) fn ungapped(cutoffs: LinkCutoffs, query_translated: bool) -> Self {
        Self {
            gap_prob: cutoffs.gap_prob,
            gap_decay_rate: BLAST_GAP_DECAY_RATE,
            cutoff_small_gap: cutoffs.cutoff_small_gap,
            cutoff_big_gap: cutoffs.cutoff_big_gap,
            query_translated,
            subject_translated: false,
        }
    }

    /// The parameters of TBLASTN's ungapped linking with one subject: `ungapped` with a protein
    /// query and a translated subject.
    pub(crate) fn ungapped_translated_subject(cutoffs: LinkCutoffs) -> Self {
        Self {
            subject_translated: true,
            ..Self::ungapped(cutoffs, false)
        }
    }

    /// The query frame of the linking a context belongs to: the context divided by `strand_factor`
    /// (link_hsps.c:513; `strand_factor` is 3 for a translated query, else 1), with the division by
    /// a constant.
    fn query_frame(&self, context: usize) -> usize {
        if self.query_translated {
            context / 3
        } else {
            context
        }
    }
}

/// What `CalculateLinkHSPCutoffs` sets for one subject of an ungapped search:
/// `link_hsp_params->cutoff_small_gap`, `cutoff_big_gap` and `gap_prob` (0 and no small gaps in a
/// search space of at most 8 windows).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LinkCutoffs {
    pub(crate) cutoff_small_gap: i32,
    pub(crate) cutoff_big_gap: i32,
    pub(crate) gap_prob: f64,
}

/// `CalculateLinkHSPCutoffs` of an ungapped search with one subject, on the parameters of
/// `BlastLinkHSPParametersNew` for an ungapped search (`gap_prob` `BLAST_GAP_PROB`,
/// `gap_decay_rate` `BLAST_GAP_DECAY_RATE`, blast_parameters.c:609-617).
/// `kbp` is the Karlin block of the context with the smallest positive Lambda
/// (`s_BlastFindSmallestLambda` on `sbp->kbp`); `avg_query_length` NCBI's average of the
/// contexts' lengths: the last context's `query_offset` plus its `query_length`, less 1, divided
/// by one more than `last_context`; `subject_length` the subject's length as the search gives it
/// (`seq_arg.seq->length`, in nucleotides for a translated subject); `db_length`
/// `BlastSeqSrcGetTotLen`; `subject_translated` `Blast_SubjectIsTranslated` (or RPS-TBLASTN);
/// `cutoff_score_min` gives `word_params->cutoff_score_min`, which NCBI reads only for a search
/// space of more than 8 windows. NCBI then multiplies both cutoffs by
/// `sbp->scale_factor` cast to `Int4`, which is 1 in the searches that call this (BLASTX and TBLASTX;
/// only an RPS search scales its matrix), so the factor is left out.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1011-1081
pub(crate) fn calculate_link_hsp_cutoffs(
    kbp: &KarlinParams,
    avg_query_length: i32,
    subject_length: i64,
    db_length: i64,
    subject_translated: bool,
    cutoff_score_min: impl FnOnce() -> i32,
) -> LinkCutoffs {
    const K_EPSILON: f64 = 1.0e-9;
    let window_size = (GAP_SIZE + OVERLAP_SIZE + 1) as i64;
    let gap_prob = BLAST_GAP_PROB;
    let gap_decay_rate = BLAST_GAP_DECAY_RATE;
    let mut query_length = avg_query_length;
    let (mut subject_length, db_length) = if subject_translated {
        ((subject_length / 3) as i32, db_length / 3)
    } else {
        (subject_length as i32, db_length)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441 (`BLAST_Nint` rounds
    // half away from zero and converts to `long`)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1036-1037 (the `(Int4)`
    // of that `long`)
    // An average query length of 0 makes the logarithm -inf: NCBI's `long` is then LONG_MIN (x86-64
    // cvttsd2si; aarch64 saturates to the same value) and its low 32 bits, the `Int4`, are 0. Rust's
    // `as i64` saturates to the same LONG_MIN and `as i32` keeps the same low bits; a direct `as i32`
    // would give i32::MIN. Finite values in range convert the same either way.
    let expected = (kbp.k * query_length as f64 * subject_length as f64).glibc_ln() / kbp.h;
    let expected_length = (expected + if expected >= 0.0 { 0.5 } else { -0.5 }) as i64 as i32;
    query_length = (query_length - expected_length).max(1);
    subject_length = (subject_length - expected_length).max(1);
    let y_variable = if db_length > subject_length as i64 {
        (db_length as f64 / subject_length as f64).glibc_ln() * kbp.k / gap_decay_rate
    } else {
        ((subject_length + expected_length) as f64 / subject_length as f64).glibc_ln() * kbp.k
            / gap_decay_rate
    };
    let search_sp = query_length as i64 * subject_length as i64;
    let mut x_variable = 0.25 * y_variable * search_sp as f64;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1061-1072
    // The `(Int4)` of `floor(log(x) / Lambda)` is NCBI's x86-64 conversion
    // (`ncbi_int4_from_double`). Where the database is not longer than the subject less the
    // expected length (one subject, for an expected length of 0), an expected length of 0 makes
    // `y_variable` and `x_variable` 0 and the logarithm -inf, and a negative one makes them
    // negative and the logarithm NaN; both convert to `INT_MIN`, so the big-gap cutoff is one more
    // than `INT_MIN` (Rust's `as i32` would give 0 for NaN, a cutoff of 1).
    let cutoff = |x: f64| ncbi_int4_from_double((x.glibc_ln() / kbp.lambda).floor()) + 1;
    let (gap_prob, cutoff_small_gap, cutoff_big_gap) = if search_sp > 8 * window_size * window_size
    {
        x_variable /= 1.0 - gap_prob + K_EPSILON;
        let cutoff_big_gap = cutoff(x_variable);
        x_variable = y_variable * (window_size * window_size) as f64;
        x_variable /= gap_prob + K_EPSILON;
        let cutoff_small_gap = cutoff_score_min().max(cutoff(x_variable));
        (gap_prob, cutoff_small_gap, cutoff_big_gap)
    } else {
        (0.0, 0, cutoff(x_variable))
    };
    LinkCutoffs {
        cutoff_small_gap,
        cutoff_big_gap,
        gap_prob,
    }
}

/// One HSP as the linking reads it (`BlastHSP`): its context, its subject frame
/// (`subject.frame`: 0 for an untranslated subject; TBLASTN's ungapped search gives the frame of
/// the translated subject, whose sign splits the query frames and orders the HSPs after the
/// context in `s_RevCompareHSPsTbn`), its query and subject ranges (0-based, the end exclusive,
/// `query.offset`/`query.end`, in the context's letters) and its raw score.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct LinkHsp {
    pub(crate) context: usize,
    pub(crate) subject_frame: i8,
    pub(crate) query_offset: i32,
    pub(crate) query_end: i32,
    pub(crate) subject_offset: i32,
    pub(crate) subject_end: i32,
    pub(crate) score: i32,
}

/// The values of a query context the linking reads: the Karlin block of the context (gapped,
/// `sbp->kbp_gap[context]`, scaled in the composition-based redo; ungapped, `sbp->kbp[context]`)
/// and the context's lengths (`query_info->contexts[context]`).
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:63-97 (LinkHSPStruct, BlastHSPLink)
/// `LinkHSPStruct`; the links are slots of `even_gap_link_hsps`'s `nodes` (slot 0 is `hp_start`, the HSP of
/// sorted position `i` is slot `i + 1`). `num` is NCBI's `Int2`: a set of more than 32,767 HSPs
/// wraps on the store of `H_hsp_num` plus 1 into `H->hsp_link.num[index]` (an `Int4` converted to
/// `Int2`, deterministic, so kept as NCBI's; the sums are `i64`, Owner 2026-10-10, as their `Int4`
/// overflow is undefined), and every later reading takes the wrapped value: the sum statistics then
/// give a set of a count below 1 the E-value 0 (`s_BlastSumP` returns 0 for `r < 1`).
#[derive(Debug, Clone, Default)]
struct Node {
    link: [Link; 2],
    num: [i16; 2],
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
    ordering_method: u8,
    evalue: f64,
    prev: Link,
    next: Link,
}

/// A slot of `even_gap_link_hsps`'s `nodes` or none (NCBI's `LinkHSPStruct*` or `NULL`), in 4 bytes: a list
/// holds far fewer than `u32::MAX` HSPs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Link(u32);

impl Link {
    const NONE: Link = Link(u32::MAX);

    fn some(slot: usize) -> Link {
        debug_assert!(slot < u32::MAX as usize);
        Link(slot as u32)
    }

    fn get(self) -> Option<usize> {
        (self != Link::NONE).then_some(self.0 as usize)
    }
}

impl Default for Link {
    fn default() -> Link {
        Link::NONE
    }
}

impl From<Option<usize>> for Link {
    fn from(slot: Option<usize>) -> Link {
        slot.map_or(Link::NONE, Link::some)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:100-109
/// `LinkHelpStruct` (`maxsum1` is read only in a test that NCBI disables with the constant 0 as
/// its condition, so it is not kept).
#[derive(Debug, Clone, Copy, Default)]
struct Helper {
    ptr: usize,
    q_off_trim: i32,
    s_off_trim: i32,
    sum: [i64; 2],
    next_larger: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:276-318
/// `s_RevCompareHSPsTbn`, with the HSPs' linking groups (`group1`, `group2`, `link_group`) in
/// place of the context and the sign test of the subject frames: `s_RevCompareHSPsTbx`
/// (link_hsps.c:331-375) for a translated query, whose `strand_factor` 3 is half of `NUM_FRAMES`.
fn rev_compare_tbn(h1: &LinkHsp, group1: usize, h2: &LinkHsp, group2: usize) -> Ordering {
    group1
        .cmp(&group2)
        .then(h2.query_offset.cmp(&h1.query_offset))
        .then(h2.query_end.cmp(&h1.query_end))
        .then(h2.subject_offset.cmp(&h1.subject_offset))
        .then(h2.subject_end.cmp(&h1.subject_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:197-224
// `s_RevCompareHSPs` (and `s_RevCompareHSPsTransl`, link_hsps.c:226-262) compares the fields that
// `s_FwdCompareHSPs` compares, in the other direction; `even_gap_link_hsps` explains why the stable
// sort by `fwd_compare` alone gives the order of the two sorts.

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:119-146
/// `s_FwdCompareHSPs`, and with a `strand_factor` of 3 `s_FwdCompareHSPsTransl`
/// (link_hsps.c:157-188).
fn fwd_compare(h1: &SortedHsp, h2: &SortedHsp) -> Ordering {
    h1.frame
        .cmp(&h2.frame)
        .then(h1.hsp.query_offset.cmp(&h2.hsp.query_offset))
        .then(h1.hsp.subject_offset.cmp(&h2.hsp.subject_offset))
}

/// The group of the HSPs that the linking links together, ascending as `s_RevCompareHSPsTbn`
/// orders them: the query frame (the context divided by `strand_factor`), then the sign of the
/// subject frame (0 for an untranslated subject; of two frames of different signs the larger frame
/// has the larger sign). Taken once for each HSP, so that the sort compares one number.
fn link_group(query_frame: usize, subject_frame: i8) -> usize {
    query_frame * 3 + (subject_frame.signum() + 1) as usize
}

/// An HSP in the order of the first sort, with what the linking reads of it on every pass: its
/// query frame (the context divided by `strand_factor`) and the terms its score adds to a set's
/// `xsum`, the score times `kbp[context]->Lambda` and `kbp[context]->logK` (the same values NCBI
/// computes on each pass, so the sums keep their bits).
#[derive(Debug, Clone, Copy)]
struct SortedHsp {
    hsp: LinkHsp,
    frame: usize,
    score_lambda: f64,
    log_k: f64,
}

/// `BLAST_LinkHsps`'s linking of a list `items` in place: each item's linking view (`hsp`) is
/// linked (`even_gap_link_hsps`), the list takes the order the linking builds, and `take` gives
/// each item its linked set's values. The caller then sorts the list by score and sets its best
/// e-value (link_hsps.c:1801-1811).
pub(crate) fn link_in_place<T>(
    items: &mut Vec<T>,
    hsp: impl Fn(&T) -> LinkHsp,
    contexts: impl Fn(usize) -> LinkContext,
    subject_length: i32,
    params: &EvenGapLinkParams,
    take: impl Fn(&mut T, &LinkedHsp),
) {
    let hsps: Vec<LinkHsp> = items.iter().map(hsp).collect();
    let mut slots: Vec<Option<T>> = std::mem::take(items).into_iter().map(Some).collect();
    items.reserve_exact(slots.len());
    even_gap_link_hsps(&hsps, &contexts, subject_length, params, &mut |linked| {
        let mut item = slots[linked.index]
            .take()
            .expect("the linking lists every HSP once");
        take(&mut item, &linked);
        items.push(item);
    });
}

/// Links the HSPs of one subject (`s_BlastEvenGapLinkHSPs` for an untranslated subject) and
/// gives `emit` the list NCBI builds, in its order: every HSP once, with its linked set's e-value
/// and size. `contexts` gives the values of a context (`LinkHsp::context`); `subject_length` is the
/// subject's length in letters (in nucleotides for a translated subject,
/// `EvenGapLinkParams::subject_translated`). `contexts` is a trait object so that one copy of the
/// linking serves every caller (a generic one was copied into each, which made the binary 35 KB
/// larger and moved the default searches' code pages).
///
/// NCBI's `qsort` calls are stable sorts here (the comparators are partial: HSPs with the same
/// start keep the order of the previous sort), as TBLASTX's port sorts them.
pub(crate) fn even_gap_link_hsps(
    hsps: &[LinkHsp],
    contexts: &dyn Fn(usize) -> LinkContext,
    subject_length: i32,
    params: &EvenGapLinkParams,
    emit: &mut dyn FnMut(LinkedHsp),
) {
    let total_number_of_hsps = hsps.len();
    if total_number_of_hsps == 0 {
        return;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:465-495
    // `link_hsp_params->gap_size` and `overlap_size` are the constants (`EvenGapLinkParams`); the
    // window of 50 is then the one every caller of `small_gap_sum_e` passes.
    let window_size = GAP_SIZE + OVERLAP_SIZE + 1;
    let trim_size = (OVERLAP_SIZE + 1) / 2;
    let gap_prob = params.gap_prob;
    let gap_decay_rate = params.gap_decay_rate;
    let cutoff = [params.cutoff_small_gap as i64, params.cutoff_big_gap as i64];
    let ignore_small_gaps = params.cutoff_small_gap == 0;
    // `link_hsp_array`: the input positions in the order of `s_RevCompareHSPsTbn` (`Tbx`).
    let frame = |h: &LinkHsp| params.query_frame(h.context);
    // Each entry holds an HSP's linking group (`link_group`, taken once for each HSP) in its high
    // half and the HSP's input position in its low half: one array of the HSPs' number, as many
    // subjects of few HSPs link one by one.
    assert!(
        u32::try_from(total_number_of_hsps).is_ok(),
        "more HSPs than 32-bit positions"
    );
    let group = |entry: u64| (entry >> 32) as usize;
    let position = |entry: u64| entry as u32 as usize;
    let mut order: Vec<u64> = hsps
        .iter()
        .enumerate()
        .map(|(i, h)| ((link_group(frame(h), h.subject_frame) as u64) << 32) | i as u64)
        .collect();
    order.sort_by(|&a, &b| {
        rev_compare_tbn(&hsps[position(a)], group(a), &hsps[position(b)], group(b))
    });
    // Slot 0 is `hp_start`; the HSP at sorted position `i` is slot `i + 1`, read from a copy in that
    // order (one load where `hsps[order[i]]` takes two).
    let sorted: Vec<SortedHsp> = order
        .iter()
        .map(|&entry| {
            let hsp = hsps[position(entry)];
            let c = contexts(hsp.context);
            SortedHsp {
                hsp,
                frame: frame(&hsp),
                score_lambda: hsp.score as f64 * c.lambda,
                log_k: c.log_k,
            }
        })
        .collect();
    let mut nodes: Vec<Node> = vec![Node::default(); total_number_of_hsps + 1];
    let hsp_of = |slot: usize| -> &LinkHsp { &sorted[slot - 1].hsp };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:510-536
    // (frame start slot, number of HSPs) of each query frame (a context, or the three contexts of a
    // strand of a translated query).
    let mut frames: Vec<(usize, usize)> = vec![(1, 0)];
    for index in 0..total_number_of_hsps {
        let slot = index + 1;
        nodes[slot].start_of_chain = false;
        frames.last_mut().expect("one frame").1 += 1;
        nodes[slot].prev = Link::from((index > 0).then_some(slot - 1));
        nodes[slot].next = Link::from((index + 1 < total_number_of_hsps).then_some(slot + 1));
        if index > 0 && group(order[index]) != group(order[index - 1]) {
            frames.last_mut().expect("one frame").1 -= 1;
            frames.push((slot, 1));
            nodes[slot - 1].next = Link::NONE;
            nodes[slot].prev = Link::NONE;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:538-552
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
    let mut lh_helper: Vec<Helper> = vec![Helper::default(); total_number_of_hsps + 5];

    for &(frame_start, frame_count) in &frames {
        link_frame(
            &mut nodes,
            &mut lh_helper,
            frame_start,
            frame_count,
            &sorted,
            contexts,
            subject_length,
            LinkConstants {
                window_size,
                trim_size,
                gap_prob,
                gap_decay_rate,
                cutoff,
                ignore_small_gaps,
                subject_translated: params.subject_translated,
            },
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:988-1001
    // The two sorts compare the same fields (the query frame, `query.offset`, `subject.offset`; the
    // offsets in the other direction), so a pair is equal under one exactly when it is equal under
    // the other. As stable sorts, the first keeps HSPs with equal fields in the slots' order and the
    // second keeps that order among them, so the second alone, on the slots' order, gives the same
    // list.
    let mut array: Vec<usize> = (1..=total_number_of_hsps).collect();
    array.sort_by(|&a, &b| fwd_compare(&sorted[a - 1], &sorted[b - 1]));

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1006-1081
    // The list NCBI builds visits, in the array's order, each HSP that is not inside a chain,
    // and after a chain's head the rest of its chain.
    let mut emitted = 0usize;
    for &slot in &array {
        let head = &nodes[slot];
        if head.linked_set && !head.start_of_chain {
            continue;
        }
        let method = usize::from(head.ordering_method);
        if head.link[method].get().is_none() {
            emit(LinkedHsp {
                index: position(order[slot - 1]),
                evalue: head.evalue,
                num: 1,
            });
            emitted += 1;
            continue;
        }
        let num_links = i32::from(head.num[method]);
        let mut current = Some(slot);
        while let Some(h) = current {
            emit(LinkedHsp {
                index: position(order[h - 1]),
                evalue: nodes[h].evalue,
                num: num_links,
            });
            emitted += 1;
            current = nodes[h].link[method].get();
        }
    }
    debug_assert_eq!(
        emitted, total_number_of_hsps,
        "the linking lists every HSP once"
    );
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
    subject_translated: bool,
}

/// The linking of one query frame (the loop of `s_BlastEvenGapLinkHSPs` over `frame_index`, from 0).
#[allow(clippy::too_many_arguments)]
fn link_frame(
    nodes: &mut [Node],
    lh_helper: &mut [Helper],
    frame_start: usize,
    frame_count: usize,
    sorted: &[SortedHsp],
    contexts: &dyn Fn(usize) -> LinkContext,
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
        subject_translated,
    } = k;
    let hsp_of = |slot: usize| -> &LinkHsp { &sorted[slot - 1].hsp };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:554-591
    const HP_START: usize = 0;
    nodes[HP_START] = Node::default();
    nodes[HP_START].next = Link::some(frame_start);
    nodes[frame_start].prev = Link::some(HP_START);
    let mut number_of_hsps = frame_count;
    // The lengths and the search space of the frame's first HSP's context, for every HSP of a
    // frame of a translated query too (NCBI's `query_context`).
    let query_context = hsp_of(frame_start).context;
    let ctx = contexts(query_context);
    let mut length_adjustment = ctx.length_adjustment;
    let query_length = (ctx.query_length - length_adjustment).max(1);
    let mut subject_length = subject_length_orig;
    if subject_translated {
        length_adjustment /= 3;
        subject_length /= 3;
    }
    let subject_length = (subject_length - length_adjustment).max(1);
    let eff_searchsp = ctx.eff_searchsp;

    lh_helper[0].ptr = HP_START;
    lh_helper[0].q_off_trim = 0;
    lh_helper[0].s_off_trim = 0;
    lh_helper[0].next_larger = 0;

    let mut first_pass = true;
    let mut path_changed = true;
    let mut h = nodes[HP_START].next.get();
    while let Some(slot) = h {
        nodes[slot].changed = true;
        h = nodes[slot].next.get();
    }

    let score_of = |slot: usize| hsp_of(slot).score as i64;
    // `H_hsp_xsum` plus the score times `kbp[H->hsp->context]->Lambda`, minus
    // `kbp[H->hsp->context]->logK`, added in C's order: the score's term first, then logK
    // subtracted.
    let add_xsum = |h_hsp_xsum: f64, slot: usize| {
        let s = &sorted[slot - 1];
        h_hsp_xsum + s.score_lambda - s.log_k
    };

    while number_of_hsps > 0 {
        let mut best: [Option<usize>; 2] = [None, None];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:603-656
        let mut use_current_max = false;
        if !first_pass {
            if !ignore_small_gaps {
                let mut max0 = -cutoff[SMALL];
                let mut max1 = -cutoff[LARGE];
                let mut h = nodes[HP_START].next.get();
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
                    h = nodes[slot].next.get();
                }
            } else {
                let mut maxscore = -cutoff[LARGE];
                let mut h = nodes[HP_START].next.get();
                while let Some(slot) = h {
                    let sum = nodes[slot].sum[LARGE];
                    if sum >= maxscore {
                        maxscore = sum;
                        best[LARGE] = Some(slot);
                    }
                    h = nodes[slot].next.get();
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
                        h = nodes[slot].link[method].get();
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:660-692
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
                h = node.next.get();
                h_index += 1;
            }

            if !ignore_small_gaps {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:693-768
                let index = SMALL;
                let mut maxscore = -cutoff[index];
                let mut h_index = 2usize;
                let mut h = nodes[HP_START].next.get();
                while let Some(slot) = h {
                    let mut h_hsp_num = 0i16;
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
                    let new_xsum = add_xsum(h_hsp_xsum, slot);
                    let new_sum = h_hsp_sum + (score_of(slot) - cutoff[index]);
                    let node = &mut nodes[slot];
                    node.sum[index] = new_sum;
                    node.num[index] = (i32::from(h_hsp_num) + 1) as i16;
                    node.link[index] = Link::from(h_hsp_link);
                    lh_helper[h_index].sum[index] = new_sum;
                    if new_sum >= maxscore {
                        maxscore = new_sum;
                        best[index] = Some(slot);
                    }
                    node.xsum[index] = new_xsum;
                    let next = node.next.get();
                    if let Some(link) = h_hsp_link {
                        nodes[link].linked_to += 1;
                    }
                    h = next;
                    h_index += 1;
                }
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:770-898
            let index = LARGE;
            let mut maxscore = -cutoff[index];
            let mut h_index = 2usize;
            let mut h = nodes[HP_START].next.get();
            while let Some(slot) = h {
                let mut h_hsp_num = 0i16;
                let mut h_hsp_sum = 0i64;
                let mut h_hsp_xsum = 0.0f64;
                let mut h_hsp_link: Option<usize> = None;
                nodes[slot].changed = true;
                let h2 = nodes[slot].link[index].get();
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
                let new_xsum = add_xsum(h_hsp_xsum, slot);
                let new_sum = h_hsp_sum + (score_of(slot) - cutoff[index]);
                let node = &mut nodes[slot];
                node.sum[index] = new_sum;
                node.num[index] = (i32::from(h_hsp_num) + 1) as i16;
                node.link[index] = Link::from(h_hsp_link);
                lh_helper[h_index].sum[index] = new_sum;
                lh_helper[h_index].next_larger = next_larger(lh_helper, h_index);
                if new_sum >= maxscore {
                    maxscore = new_sum;
                    best[index] = Some(slot);
                }
                node.xsum[index] = new_xsum;
                let next = node.next.get();
                if let Some(link) = h_hsp_link {
                    nodes[link].linked_to += 1;
                }
                h = next;
                h_index += 1;
            }
            path_changed = false;
            first_pass = false;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:891-952
        let large_gap_e = |nodes: &[Node], best: usize| {
            let num = nodes[best].num[LARGE];
            ncbi_large_gap_sum_e(
                num,
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
            nodes[b0].sum[SMALL] += i64::from(nodes[b0].num[SMALL]) * cutoff[SMALL];
            let num0 = nodes[b0].num[SMALL];
            prob[SMALL] = small_gap_sum_e(
                window_size,
                num0,
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
            nodes[b1].sum[LARGE] += i64::from(nodes[b1].num[LARGE]) * cutoff[LARGE];
            prob[LARGE] = large_gap_e(nodes, b1);
            ordering_method = LARGE;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:948-980
        let chosen = best[ordering_method].expect("the chosen ordering has a best set");
        nodes[chosen].start_of_chain = true;
        nodes[chosen].evalue = prob[ordering_method];
        let linked_set = nodes[chosen].link[ordering_method].get().is_some();
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
            node.ordering_method = ordering_method as u8;
            node.evalue = prob[ordering_method];
            let (prev, next) = (node.prev, node.next);
            if let Some(n) = next.get() {
                nodes[n].prev = prev;
            }
            if let Some(p) = prev.get() {
                nodes[p].next = next;
            }
            number_of_hsps -= 1;
            h = nodes[slot].link[ordering_method].get();
        }
    }
}

/// `next_larger` of `lh_helper[h_index]`: the closest earlier entry whose large-gap sum is
/// larger.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:676-688 (the same at 864-873)
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
#[path = "link_hsps_tests.rs"]
mod tests;
