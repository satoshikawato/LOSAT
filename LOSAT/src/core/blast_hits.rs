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
/// TBLASTX's ungapped re-evaluation (`blast_hits.c:2700-2705`) refuse there
/// (`x_untranslated_frame_refusal`).
#[inline]
pub(crate) fn target_translation_never_made(subject_length: usize) -> bool {
    subject_length / 3 <= 3
}

/// The refusal of `program`'s `step` (the NCBI step that reads the translation) of an alignment
/// on a subject whose frames NCBI never translates (`target_translation_never_made`): NCBI BLAST+
/// 2.17.0 reads from address 1 on (NULL + 1) and dies with SIGSEGV (PD-LOSAT-NCBI-DEFECTS 1.46 for
/// TBLASTN, 1.51 for TBLASTX, rule 3).
#[cold]
#[inline(never)]
pub(crate) fn x_untranslated_frame_refusal(program: &str, step: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "{program}: {step} of an alignment on a subject (or its -subject_loc range) of fewer than 12 bases, where NCBI BLAST+ reads the translation of the subject's frame, which it never makes for a subject this short, and crashes; LOSAT does not support it"
    )
}

#[cfg(test)]
#[path = "blast_hits_tests.rs"]
mod tests;
