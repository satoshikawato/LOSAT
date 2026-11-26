//! HSP Storage and Management
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c
//!
//! The filters of `s_FilterBlastResults` that every program shares (`-max_hsps`,
//! `-qcov_hsp_perc`), `Blast_HSPTest` (NCBI's `Int4` counts, widened to `i64`), and the frames
//! that `Blast_HSPGetTargetTranslation` never translates, with their refusal (TBLASTN, TBLASTX).

/// Whether an HSP fails NCBI's percent identity or minimal hit length test (`Blast_HSPTest`,
/// `s_HSPTest`: the HSP is deleted), for every program. The counts are NCBI's `Int4`s widened to
/// `i64`, so that the programs that count in `usize` pass their lengths unchanged (every `Int4`
/// is the same `i64` and the same `double`). The alignment length is negative for an ungapped
/// HSP whose end lies before its start (BLASTX's one-hit word finder's HSPs of length below 0
/// with `-window_size 0`, `aa_ungapped.c:1054,1083`): `align_length < min_hit_length` holds
/// there for every `min_hit_length` of 0 or more. The identity test multiplies the length by
/// `percent_identity` in `double`, as NCBI does. Neither test needs a guard for its disabled
/// value: with `-perc_identity` 0 (NCBI's range is 0 to 100, `blast_args.cpp:188-192`) and the
/// `min_hit_length` of 0 that NCBI's applications always leave (`calloc`, no option), both are
/// false for every length of 0 or more.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1446,1448
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-1001,1027-1032
pub fn hsp_test(
    num_ident: i64,
    align_length: i64,
    percent_identity: f64,
    min_hit_length: i64,
) -> bool {
    // The counts are lengths of sequences in memory: far below 2^53, exact in `f64`.
    num_ident as f64 * 100.0 < align_length as f64 * percent_identity
        || align_length < min_hit_length
}

/// The query coverage of an HSP in percent: its query extent (`offset` to `end`, end
/// exclusive, in the letters of its context) over the length of its context, plus 0.5 below
/// 99 (`Blast_HSPGetQueryCoverage`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1034-1044
pub fn hsp_query_coverage(offset: i64, end: i64, query_length: i64) -> f64 {
    let mut pct = 0.0;
    if query_length > 0 {
        pct = 100.0 * (end - offset) as f64 / query_length as f64;
        if pct < 99.0 {
            pct += 0.5;
        }
    }
    pct
}

/// Whether an HSP covers less than `min_pct` percent of its query context
/// (`Blast_HSPQueryCoverageTest`: the HSP is deleted).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1045-1051
pub fn hsp_query_coverage_test(min_pct: f64, offset: i64, end: i64, query_length: i64) -> bool {
    hsp_query_coverage(offset, end, query_length) < min_pct
}

/// `Blast_HSPListReapByQueryCoverage` over the HSPs of one list, in their order: the HSPs
/// whose query coverage (`hsp_query_coverage_test` of the `extent` of each: query offset,
/// query end and the length of its context) is below `min_pct` are removed. A `min_pct` of
/// 0 is NCBI's unset option and removes nothing. Returns whether an HSP was removed (the
/// caller then sets the list's best e-value again, and frees a list left empty).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2010-2045
pub fn reap_by_query_coverage<T>(
    hsps: &mut Vec<T>,
    min_pct: f64,
    extent: impl Fn(&T) -> (i64, i64, i64),
) -> bool {
    if hsps.is_empty() || min_pct == 0.0 {
        return false;
    }
    let count = hsps.len();
    hsps.retain(|hsp| {
        let (offset, end, query_length) = extent(hsp);
        !hsp_query_coverage_test(min_pct, offset, end, query_length)
    });
    hsps.len() != count
}

/// `Blast_TrimHSPListByMaxHsps`: the first `max_hsps` HSPs of a list, in its order, are
/// kept (0 keeps all). The list's best e-value is not set again.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2049-2069
pub fn trim_by_max_hsps<T>(hsps: &mut Vec<T>, max_hsps: usize) {
    if max_hsps == 0 || hsps.len() <= max_hsps {
        return;
    }
    hsps.truncate(max_hsps);
}

/// The query range around a kept HSP within which `-subject_besthit` removes the later HSPs of
/// its context, in the letters of the context: 3 for every program (NCBI's protein and
/// nucleotide defaults are both 3).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:1211-1212
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1962-1974
pub const SUBJECT_BESTHIT_MAX_RANGE_DIFF: i64 = 3;

/// The fields of an HSP that `-subject_besthit` reads: its context (NCBI's context index: the
/// query, its strand for BLASTN and its frame for BLASTX and TBLASTX), the query offset and end
/// in the letters of that context (NCBI's `Int4`s, widened; an HSP can lie outside its context,
/// `blastn::hsp::BlastnHsp::outside_query_context`, or end before its start, BLASTX's and
/// TBLASTX's one-hit HSPs with `-window_size 0`), and, for BLASTN's pass over the other strand,
/// the sign of the query frame and the length of the context.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BestHitExtent {
    pub context: i64,
    pub query_frame: i32,
    pub offset: i64,
    pub end: i64,
    pub query_length: i64,
}

/// `Blast_HSPListSubjectBestHit` over the HSPs of one subject's list, in their order (NCBI's
/// list is sorted by score where it is called; the filter reads only the order). Each kept HSP,
/// from the first on, removes every later HSP of its context whose query range lies within its
/// own widened by `SUBJECT_BESTHIT_MAX_RANGE_DIFF` on both sides (the start not below 0; an end
/// that would be negative is the HSP's own end). For BLASTN (`blastn`) a second pass over the
/// kept HSPs then removes, from the other strand's context of the same query, the later HSPs
/// whose query range lies within the kept HSP's range reflected onto that strand and widened the
/// same way (no bound at 0 there). The kept HSPs keep their order. The list's best e-value is
/// not set again (NCBI does not), and an empty list stays empty.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2537-2606
pub fn subject_besthit<T>(hsps: &mut Vec<T>, blastn: bool, extent: impl Fn(&T) -> BestHitExtent) {
    if hsps.len() <= 1 {
        return;
    }
    let range_diff = SUBJECT_BESTHIT_MAX_RANGE_DIFF;
    let mut extents: Vec<BestHitExtent> = hsps.iter().map(&extent).collect();
    let mut removed = vec![false; extents.len()];
    for i in 0..extents.len() {
        if removed[i] {
            continue;
        }
        let kept = extents[i];
        let low = (kept.offset - range_diff).max(0);
        let mut high = kept.end + range_diff;
        if high < 0 {
            high = kept.end;
        }
        for (j, later) in extents.iter().enumerate().skip(i + 1) {
            if !removed[j]
                && later.context == kept.context
                && later.offset >= low
                && later.end <= high
            {
                removed[j] = true;
            }
        }
    }
    retain_unremoved(hsps, &removed);
    if !blastn || hsps.len() <= 1 {
        return;
    }
    extents.clear();
    extents.extend(hsps.iter().map(&extent));
    removed.clear();
    removed.resize(extents.len(), false);
    for i in 0..extents.len() {
        if removed[i] {
            continue;
        }
        let kept = extents[i];
        let other_strand = if kept.query_frame > 0 {
            kept.context + 1
        } else {
            kept.context - 1
        };
        let high = kept.query_length - (kept.offset - range_diff);
        let low = kept.query_length - (kept.end + range_diff);
        for (j, later) in extents.iter().enumerate().skip(i + 1) {
            if !removed[j]
                && later.context == other_strand
                && later.offset >= low
                && later.end <= high
            {
                removed[j] = true;
            }
        }
    }
    retain_unremoved(hsps, &removed);
}

/// Keeps the items of `items` whose mark in `removed` is false, in their order.
fn retain_unremoved<T>(items: &mut Vec<T>, removed: &[bool]) {
    let mut marks = removed.iter();
    items.retain(|_| !marks.next().copied().unwrap_or(false));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1284;
// c++/src/algo/blast/core/blast_hits.c:1162-1164,1234-1235
/// Whether NCBI's `Blast_HSPGetTargetTranslation` never translates a frame of a subject of
/// `subject_length` bases (the searched range, `-subject_loc`: NCBI's `subject_blk->length`): its
/// test is false at the first call (`start` and `stop` 0) for fewer than 12 bases and stays so,
/// and it returns the translation it never made (`calloc`'d to NULL) plus 1, from which every
/// caller reads next. TBLASTN's traceback reads (`tblastn::search_gapped::TargetTranslation`) and
/// the ungapped re-evaluation of TBLASTX and TBLASTN `-ungapped` (`blast_hits.c:2700-2705`) refuse
/// there (`x_untranslated_frame_refusal`).
#[inline]
pub(crate) fn target_translation_never_made(subject_length: usize) -> bool {
    subject_length / 3 <= 3
}

/// The refusal of `program`'s `step` (the NCBI step that reads the translation) of an alignment
/// on a subject whose frames NCBI never translates (`target_translation_never_made`): NCBI BLAST+
/// 2.17.0 reads from address 1 on (NULL + 1) and dies with SIGSEGV (PD-LOSAT-NCBI-DEFECTS 1.46 for
/// TBLASTN, 1.79 for TBLASTN `-ungapped`, 1.51 for TBLASTX, rule 3).
#[cold]
#[inline(never)]
pub(crate) fn x_untranslated_frame_refusal(program: &str, step: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "{program}: {step} of an alignment on a subject (or its -subject_loc range) of fewer than 12 bases, where NCBI BLAST+ reads the translation of the subject's frame, which it never makes for a subject this short, and crashes; LOSAT does not support it"
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1533-1537
/// Maximum diagonal separation between an earlier HSP's end and an incoming HSP's start for
/// `Blast_HSPListsMerge` to attempt a merge across a subject-chunk overlap (`OVERLAP_DIAG_CLOSE`).
const OVERLAP_DIAG_CLOSE: i64 = 10;

/// Returned by [`merge_two_chunk_hsps`] and [`merge_subject_chunk_hsps`] when two HSPs in a chunk
/// overlap have query lengths summing to 0 or a density-scaled merged score outside NCBI's `Int4`,
/// where NCBI's `s_BlastMergeTwoHSPs` divides by zero or converts an out-of-range float to `Int4`
/// (undefined behavior: zero length sum is reachable with `-ungapped -window_size 0` one-hit HSPs
/// of length 0 or less; out-of-range scaled score is reachable when merging very long high-scoring
/// alignments whose summed score exceeds `Int4::MAX`, e.g. in BLASTN).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ChunkMergeUndefined;

/// The fields of a preliminary HSP read and updated by NCBI's split-subject
/// `Blast_HSPListsMerge` and `s_BlastMergeTwoHSPs`.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-143
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1464-1537,2857-3035
pub(crate) trait ChunkMergeHsp {
    fn merge_context(&self) -> i32;
    fn merge_subject_frame(&self) -> i32;
    fn merge_q_start(&self) -> i32;
    fn merge_q_end(&self) -> i32;
    fn merge_s_start(&self) -> i32;
    fn merge_s_end(&self) -> i32;
    fn merge_score(&self) -> i32;
    fn merge_q_gapped_start(&self) -> i32;
    fn merge_s_gapped_start(&self) -> i32;
    fn set_merge_bounds(&mut self, q_start: i32, q_end: i32, s_start: i32, s_end: i32);
    fn set_merge_gapped_start(&mut self, q_gapped_start: i32, s_gapped_start: i32);
    fn set_merge_score(&mut self, score: i32);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1464-1537
// `s_BlastMergeTwoHSPs`: for an ungapped merge requires equal diagonals; always requires equal
// subject frames and that either endpoint of `new` lies within the closed rectangle of `old`.
// Expands `old`'s endpoints, adopts `new`'s gapped start if `new` has strictly higher score, and
// sets `old`'s score to the larger of its score and the truncated density-scaled query span.
pub(crate) fn merge_two_chunk_hsps<H: ChunkMergeHsp>(
    old: &mut H,
    new: &H,
    allow_gap: bool,
) -> Result<bool, ChunkMergeUndefined> {
    let (q_off1, q_end1) = (old.merge_q_start(), old.merge_q_end());
    let (s_off1, s_end1) = (old.merge_s_start(), old.merge_s_end());
    let (q_off2, q_end2) = (new.merge_q_start(), new.merge_q_end());
    let (s_off2, s_end2) = (new.merge_s_start(), new.merge_s_end());
    if !allow_gap
        && i64::from(s_off1) - i64::from(s_off2) - i64::from(q_off1) + i64::from(q_off2) != 0
    {
        return Ok(false);
    }
    if old.merge_subject_frame() != new.merge_subject_frame() {
        return Ok(false);
    }
    let inside = |q: i32, s: i32| q >= q_off1 && q <= q_end1 && s >= s_off1 && s <= s_end1;
    if !inside(q_off2, s_off2) && !inside(q_end2, s_end2) {
        return Ok(false);
    }
    let length_sum =
        (i64::from(q_end1) - i64::from(q_off1)) + (i64::from(q_end2) - i64::from(q_off2));
    if length_sum == 0 {
        return Err(ChunkMergeUndefined);
    }
    let score_sum = i64::from(old.merge_score()) + i64::from(new.merge_score());
    let score_density = (score_sum as f64) / (length_sum as f64);
    let q_start = q_off1.min(q_off2);
    let s_start = s_off1.min(s_off2);
    let q_end = q_end1.max(q_end2);
    let s_end = s_end1.max(s_end2);
    let scaled = score_density * ((i64::from(q_end) - i64::from(q_start)) as f64);
    let truncated = scaled.trunc();
    if !truncated.is_finite() || truncated < f64::from(i32::MIN) || truncated > f64::from(i32::MAX)
    {
        return Err(ChunkMergeUndefined);
    }
    old.set_merge_bounds(q_start, q_end, s_start, s_end);
    if new.merge_score() > old.merge_score() {
        old.set_merge_gapped_start(new.merge_q_gapped_start(), new.merge_s_gapped_start());
        old.set_merge_score(new.merge_score());
    }
    old.set_merge_score((scaled as i32).max(old.merge_score()));
    Ok(true)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2749-2807,2857-3035
// Split-subject `Blast_HSPListsMerge` (`contexts_per_query` < 0) followed by
// `s_BlastHSPListsCombineByScore`: partitions `combined` by `s_end > split_offset` and `incoming`
// by `s_start < split_offset + overlap` using NCBI's swap displacement order, merges compatible
// pairs whose end/start diagonals differ by less than `OVERLAP_DIAG_CLOSE`, removes consumed
// incoming HSPs preserving order, and combines the two lists in `ScoreCompareHSPs` order (stable
// sort when out of order, matching glibc's `qsort` in pinned NCBI BLAST+ 2.17.0).
/// Merges `incoming` (already offset-adjusted by `Blast_HSPListAdjustOffsets`) into `combined`
/// across a subject chunk split at `split_offset` with overlap `overlap` (`Blast_HSPListsMerge`),
/// draining `incoming` so callers can reuse its allocation across chunks.
pub(crate) fn merge_subject_chunk_hsps<H: ChunkMergeHsp>(
    combined: &mut Vec<H>,
    incoming: &mut Vec<H>,
    split_offset: i32,
    overlap: i32,
    hsp_num_max: usize,
    allow_gap: bool,
    compare_score: impl Fn(&H, &H) -> std::cmp::Ordering,
) -> Result<(), ChunkMergeUndefined> {
    merge_subject_chunk_hsps_with_sort(
        combined,
        incoming,
        split_offset,
        overlap,
        hsp_num_max,
        allow_gap,
        &compare_score,
        |hsps| {
            if hsps
                .windows(2)
                .any(|pair| compare_score(&pair[0], &pair[1]).is_gt())
            {
                hsps.sort_by(&compare_score);
            }
        },
    )
}

/// [`merge_subject_chunk_hsps`] with a caller-supplied stable score sorter (`sort_by_score`,
/// used by TBLASTX's in-place index-replay sort).
#[allow(clippy::too_many_arguments)]
pub(crate) fn merge_subject_chunk_hsps_with_sort<H: ChunkMergeHsp>(
    combined: &mut Vec<H>,
    incoming: &mut Vec<H>,
    split_offset: i32,
    overlap: i32,
    hsp_num_max: usize,
    allow_gap: bool,
    compare: impl Fn(&H, &H) -> std::cmp::Ordering,
    sort_by_score: impl Fn(&mut [H]),
) -> Result<(), ChunkMergeUndefined> {
    if incoming.is_empty() {
        return Ok(());
    }
    if combined.is_empty() {
        std::mem::swap(combined, incoming);
        return Ok(());
    }
    let split = i64::from(split_offset);
    let strip_end = split + i64::from(overlap);

    let mut hspcnt1 = 0usize;
    for index1 in 0..combined.len() {
        if i64::from(combined[index1].merge_s_end()) > split {
            combined.swap(hspcnt1, index1);
            hspcnt1 += 1;
        }
    }
    let mut hspcnt2 = 0usize;
    for index2 in 0..incoming.len() {
        if i64::from(incoming[index2].merge_s_start()) < strip_end {
            incoming.swap(hspcnt2, index2);
            hspcnt2 += 1;
        }
    }

    let mut result = Ok(());
    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut deleted = vec![false; hspcnt2];
        'merge: for index1 in 0..hspcnt1 {
            for index2 in 0..hspcnt2 {
                if deleted[index2]
                    || incoming[index2].merge_context() != combined[index1].merge_context()
                {
                    continue;
                }
                let end_diag = i64::from(combined[index1].merge_q_end())
                    - i64::from(combined[index1].merge_s_end());
                let start_diag = i64::from(incoming[index2].merge_q_start())
                    - i64::from(incoming[index2].merge_s_start());
                if (end_diag - start_diag).abs() < OVERLAP_DIAG_CLOSE {
                    match merge_two_chunk_hsps(&mut combined[index1], &incoming[index2], allow_gap)
                    {
                        Ok(true) => deleted[index2] = true,
                        Ok(false) => {}
                        Err(undefined) => {
                            result = Err(undefined);
                            break 'merge;
                        }
                    }
                }
            }
        }
        let mut index2 = 0usize;
        incoming.retain(|_| {
            let keep = index2 >= hspcnt2 || !deleted[index2];
            index2 += 1;
            keep
        });
    }

    let retained = combined
        .len()
        .saturating_add(incoming.len())
        .min(hsp_num_max);
    if retained == combined.len() + incoming.len() {
        combined.append(incoming);
        sort_by_score(combined);
        return result;
    }
    sort_by_score(combined);
    sort_by_score(incoming);
    let mut left = combined.drain(..).peekable();
    let mut right = incoming.drain(..).peekable();
    let mut selected = Vec::with_capacity(retained);
    while selected.len() < retained {
        let take_left = match (left.peek(), right.peek()) {
            (Some(a), Some(b)) => compare(a, b).is_le(),
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => break,
        };
        if take_left {
            selected.push(left.next().expect("peeked combined HSP"));
        } else {
            selected.push(right.next().expect("peeked incoming HSP"));
        }
    }
    drop(left);
    drop(right);
    *combined = selected;
    result
}

#[cfg(test)]
#[path = "blast_hits_tests.rs"]
mod tests;
