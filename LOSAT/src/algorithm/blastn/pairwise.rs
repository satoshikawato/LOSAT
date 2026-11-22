//! The BLASTN pairwise report (outfmt 0): the displayed alignment of every HSP of the final
//! hit list, and the per-query statistics that the report prints.

use std::iter::repeat_n;

use crate::blastinput::fasta_reader::FastaRecord;
use anyhow::{ensure, Context, Result};

use crate::common::GapEditOp;
#[cfg(test)]
use crate::core::blast_encoding::encode_iupac_to_blastna;
use crate::report::pairwise::PairwiseHit;
use crate::stats::{
    compute_karlin_params_ungapped, score_freq_profile_from_probabilities, KarlinParams,
};
use crate::utils::dust::MaskedInterval;

use super::alignment::build_blastna_matrix;
use super::hsp::{BlastnHitList, BlastnHsp, BlastnHspList};
use super::lookup::{reverse_complement, reverse_complement_in_place};
use super::x_subject_store::SubjectLetters;

/// The regions that the report prints in lowercase.
pub(crate) struct DisplayMasks<'a> {
    /// Per query, on its plus strand of the searched letters: DUST and, with
    /// `-lcase_masking`, the input lowercase.
    pub query: &'a [Vec<MaskedInterval>],
    /// Where each query's searched letters lie in its record (`-query_loc`).
    pub query_offsets: &'a crate::blastinput::seq_range::Placements,
    /// Per subject, with `-lcase_masking`: the input lowercase of the subject's record.
    pub subject: Option<&'a [Vec<MaskedInterval>]>,
    /// Per reported database sequence (by rank), with `-db_soft_mask`: its stored masks of the
    /// search's algorithm in the order they are stored (`BlastnDatabase::display_masks`).
    pub database: Option<&'a [Vec<MaskedInterval>]>,
}

impl<'a> Rows<'a> {
    /// The rows' inputs of a report call of `queries` and `subjects`, its shown masks made
    /// once for every window of its queries (`pairwise_hits`).
    // Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
    // shared post-processing makes every BLASTN run compile it in Wasm hosts.
    #[inline(never)]
    pub(crate) fn new(
        queries: &'a [FastaRecord],
        subjects: SubjectLetters<'a>,
        masks: &DisplayMasks<'_>,
    ) -> Self {
        let query_masks: Vec<Vec<MaskedInterval>> = masks
            .query
            .iter()
            .zip(queries)
            .enumerate()
            .map(|(q_idx, (query_masks, query))| {
                shown_query_masks(
                    query_masks,
                    query.seq().len(),
                    masks.query_offsets.offset(q_idx),
                )
            })
            .collect();
        // A subject set's lowercase is merged first; a database's masks are kept as stored
        // (`shown_subject_masks` merges what it keeps).
        let subject_masks: Option<Vec<Vec<MaskedInterval>>> = match (masks.database, masks.subject)
        {
            (Some(database), _) => Some(database.to_vec()),
            (None, Some(subject)) => Some(subject.iter().map(|masks| merged(masks)).collect()),
            (None, None) => None,
        };
        Rows {
            queries,
            subjects,
            query_masks,
            subject_masks,
        }
    }
}

/// Every HSP of the final hit lists `hit_lists` (those of all queries of `rows`, or of a window
/// of them), in the order of the tabular formats (queries, then the subjects and HSPs of each
/// query as ranked), with its displayed rows.
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
pub(crate) fn pairwise_hits(
    hit_lists: &[Option<BlastnHitList>],
    rows: &Rows<'_>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<Vec<PairwiseHit>> {
    let hsp_lists = hit_lists
        .iter()
        .flatten()
        .flat_map(|hit_list| &hit_list.hsplist_array);
    // EXPERIMENT (LOSAT_X_PAIRWISEPAR): each subject's list of each query on the pool, the rows
    // collected in the list order (`report::x_pairwise_par`); the same function of each list.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if let Some(pool) = x_pool {
        use rayon::prelude::*;
        let lists: Vec<&BlastnHspList> = hsp_lists.collect();
        let rows_of_lists: Vec<Result<Vec<PairwiseHit>>> = pool.install(|| {
            lists
                .par_iter()
                .map(|hsp_list| {
                    let mut hits = Vec::with_capacity(hsp_list.hsps.len());
                    rows.list(hsp_list, &mut hits)?;
                    Ok(hits)
                })
                .collect()
        });
        let mut hits = Vec::with_capacity(lists.iter().map(|list| list.hsps.len()).sum());
        for list_rows in rows_of_lists {
            hits.append(&mut list_rows?);
        }
        return Ok(hits);
    }
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let _ = x_pool;
    let mut hits = Vec::new();
    for hsp_list in hsp_lists {
        rows.list(hsp_list, &mut hits)?;
    }
    Ok(hits)
}

/// The windows of the queries `0..num_queries` whose outfmt 0 rows are made and written
/// together (`x_pairwise_par::query_windows` of the estimates of their HSPs' lengths in a report
/// of `line_length` columns a line, the report's `PairwiseConfig::line_length`, before any row is
/// made, with the budget of `x_pool`, `x_pairwise_par::window_bytes`).
#[inline(never)]
pub(crate) fn pairwise_windows(
    hit_lists: &[Option<BlastnHitList>],
    num_queries: usize,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    line_length: usize,
) -> Vec<std::ops::Range<usize>> {
    use crate::report::x_pairwise_par::{query_estimate_of, query_windows, window_bytes};
    let estimates: Vec<usize> = (0..num_queries)
        .map(|q_idx| {
            query_estimate_of(
                hit_lists
                    .get(q_idx)
                    .and_then(Option::as_ref)
                    .into_iter()
                    .flat_map(|hit_list| &hit_list.hsplist_array)
                    .flat_map(|hsp_list| &hsp_list.hsps)
                    .map(|hsp| hsp.length),
                line_length,
            )
        })
        .collect();
    query_windows(&estimates, window_bytes(x_pool))
}

/// The displayed rows of one HSP (the query on its plus strand, the subject reverse-complemented
/// with it for a minus-strand HSP, upper case, no masking) for BLAST XML, which recounts its
/// identities from them (`algorithm::blastn::blast_xml`). The same function of the HSP as the
/// outfmt 0 rows, without a `PairwiseHit`, written into two buffers that the caller reuses.
pub(crate) fn xml_rows_into(
    hsp: &BlastnHsp,
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_row: &mut Vec<u8>,
    subject_row: &mut Vec<u8>,
) -> Result<()> {
    let query = queries
        .get(hsp.q_idx as usize)
        .context("BLASTN HSP of an unknown query")?;
    let s_idx = hsp.s_idx as usize;
    subjects
        .records()
        .get(s_idx)
        .context("BLASTN HSP of an unknown subject")?;
    displayed_rows_into(
        hsp,
        query.seq(),
        subjects.len(s_idx),
        |from, to| subjects.upper_letters(s_idx, from, to),
        query_row,
        subject_row,
    )
}

/// What the displayed rows of an HSP read besides the HSP (`pairwise_hits`).
pub(crate) struct Rows<'a> {
    queries: &'a [FastaRecord],
    subjects: SubjectLetters<'a>,
    /// Per query, the shown masks (`shown_query_masks`).
    query_masks: Vec<Vec<MaskedInterval>>,
    /// Per subject, the masks that `shown_subject_masks` shows with its HSPs: of a subject set's
    /// record with `-lcase_masking`, its merged lowercase (`DisplayMasks::subject`); of a reported
    /// database sequence with `-db_soft_mask`, its stored masks as stored
    /// (`DisplayMasks::database`).
    subject_masks: Option<Vec<Vec<MaskedInterval>>>,
}

impl Rows<'_> {
    /// The shown masks of query `q_idx` (`shown_query_masks`), the query's `TMaskedQueryRegions`
    /// in record positions (the query-anchored reports, `algorithm::blastn::anchored`).
    pub(crate) fn query_masks(&self, q_idx: usize) -> &[MaskedInterval] {
        self.query_masks.get(q_idx).map_or(&[], Vec::as_slice)
    }

    /// The masks of the subject of the HSPs `hsps` (one subject's list of a query) that NCBI's
    /// `GetSubjectMasks` gives with them (`subject_mask_regions`), in their order and unmerged:
    /// the query-anchored reports place each one on the first row that meets it
    /// (`query_anchored::fill_loc_list`). None without subject masks.
    pub(crate) fn subject_masks(&self, hsps: &[BlastnHsp]) -> Vec<MaskedInterval> {
        match (self.subject_masks.as_deref(), hsps.first()) {
            (Some(subject_masks), Some(hsp)) => {
                subject_mask_regions(&subject_masks[hsp.s_idx as usize], hsps).collect()
            }
            _ => Vec::new(),
        }
    }

    /// The rows of every HSP of one subject's list, appended to `hits` in the list's order.
    fn list(&self, hsp_list: &BlastnHspList, hits: &mut Vec<PairwiseHit>) -> Result<()> {
        let (queries, subjects) = (self.queries, self.subjects);
        let list_subject_masks = self.subject_masks.as_deref().and_then(|subject_masks| {
            let s_idx = hsp_list.hsps.first()?.s_idx as usize;
            Some(shown_subject_masks(&subject_masks[s_idx], &hsp_list.hsps))
        });
        for hsp in &hsp_list.hsps {
            let query = queries
                .get(hsp.q_idx as usize)
                .context("BLASTN HSP of an unknown query")?;
            let s_idx = hsp.s_idx as usize;
            let subject = subjects
                .records()
                .get(s_idx)
                .context("BLASTN HSP of an unknown subject")?;
            let subject_length = subjects.len(s_idx);
            let (mut query_row, mut subject_row) =
                displayed_rows(hsp, query.seq(), subject_length, |from, to| {
                    subjects.upper_letters(s_idx, from, to)
                })?;
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1794-1798
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2127-2140
            // The identities that outfmt 0 shows are those of the displayed rows (the stored
            // letters in upper case, before the masks' lower case), counted again (series NP,
            // M4a-2b decision H6): with `-db_hard_mask` the engine's `num_ident` counts a masked
            // letter ('B' in the traceback) as a mismatch. Elsewhere the two counts are equal.
            let displayed_ident = query_row
                .iter()
                .zip(&subject_row)
                .filter(|(query, subject)| query == subject)
                .count();
            lowercase_masked(
                &mut query_row,
                hsp.q_start,
                1,
                &self.query_masks[hsp.q_idx as usize],
            );
            if let Some(subject_masks) = &list_subject_masks {
                let step = if hsp.s_start > hsp.s_end { -1 } else { 1 };
                lowercase_masked(&mut subject_row, hsp.s_start, step, subject_masks);
            }
            let gaps = query_row
                .iter()
                .chain(&subject_row)
                .filter(|&&residue| residue == b'-')
                .count();
            let mut hit = hsp.clone().into_hit();
            hit.num_ident = displayed_ident;
            hits.push(PairwiseHit {
                hit,
                query_seq: Some(String::from_utf8(query_row).context("BLASTN query row")?),
                subject_seq: Some(String::from_utf8(subject_row).context("BLASTN subject row")?),
                query_frame: None,
                subject_frame: None,
                positives: None,
                gaps: Some(gaps),
                subject_length: Some(subject_length),
                // The title after its first word; the reports take the subject's title
                // bytes by subject index.
                subject_title: subject
                    .title
                    .iter()
                    .position(|&byte| byte == b' ')
                    .map(|space| String::from_utf8_lossy(&subject.title[space + 1..]).into_owned()),
                comp_adjust_method: None,
                sum_n: None,
                subject_stored_gaps: Vec::new(),
                query_stored_gaps: Vec::new(),
                query_case: Default::default(),
                // The Seq-align's `num_ident` (`seqalign::query_records` writes the same
                // `hsp.num_ident`): the engine's count, which "Max ident" and the sorts read.
                seqalign_num_ident: hsp.num_ident,
            });
        }
        Ok(())
    }
}

// The rows of one HSP as NCBI shows them: the query on its plus strand.
//
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1852-1858
// NCBI reference (598d8ae6): c++/include/algo/blast/core/gapinfo.h:43-61
// The HSP offsets and its edit script are those of the query context (the reverse
// complement of the query for a minus-strand HSP) against the plus-strand subject.
// Reversing the alignment of a minus-strand HSP reverse-complements both rows. The
// residues are the input letters in uppercase (NCBI's IUPAC sequence data); masking is
// applied afterwards. Only the aligned segments are read: the segment of the minus-strand
// context `[qs, qe)` is the reverse complement of query `[len - qe, len - qs)`. The subject
// has `subject_length` letters, and `subject_letters(from, to)` gives its letters
// `[from, to)` in uppercase.
fn displayed_rows(
    hsp: &BlastnHsp,
    query: &[u8],
    subject_length: usize,
    subject_letters: impl FnOnce(usize, usize) -> Vec<u8>,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let (mut query_row, mut subject_row) = (Vec::new(), Vec::new());
    displayed_rows_into(
        hsp,
        query,
        subject_length,
        subject_letters,
        &mut query_row,
        &mut subject_row,
    )?;
    Ok((query_row, subject_row))
}

/// `displayed_rows` into `query_row` and `subject_row` (cleared first; their capacity is kept),
/// reversed in place for a minus-strand HSP.
fn displayed_rows_into(
    hsp: &BlastnHsp,
    query: &[u8],
    subject_length: usize,
    subject_letters: impl FnOnce(usize, usize) -> Vec<u8>,
    query_row: &mut Vec<u8>,
    subject_row: &mut Vec<u8>,
) -> Result<()> {
    context_rows_into(
        hsp,
        query,
        subject_length,
        subject_letters,
        query_row,
        subject_row,
    )?;
    if hsp.query_frame < 0 {
        reverse_complement_in_place(query_row);
        reverse_complement_in_place(subject_row);
    }
    Ok(())
}

/// The two rows in the direction of the HSP's query context (cleared first; their capacity is
/// kept): a minus-strand HSP's rows are still to be reverse-complemented.
fn context_rows_into(
    hsp: &BlastnHsp,
    query: &[u8],
    subject_length: usize,
    subject_letters: impl FnOnce(usize, usize) -> Vec<u8>,
    query_row: &mut Vec<u8>,
    subject_row: &mut Vec<u8>,
) -> Result<()> {
    query_row.clear();
    subject_row.clear();
    // An HSP outside its query context (a negative `Int4` offset) never reaches a report
    // (`run::hsp_outside_query_context`). The end is not below the offset, so neither cast
    // wraps; the range test below would refuse one that did.
    ensure!(
        hsp.internal_q_offset_0 >= 0,
        "BLASTN HSP outside its sequences"
    );
    let (qs, qe) = (
        hsp.internal_q_offset_0 as usize,
        hsp.internal_q_end_0 as usize,
    );
    let (ss, se) = (hsp.internal_s_offset_0, hsp.internal_s_end_0);
    ensure!(
        qs <= qe && qe <= query.len() && ss <= se && se <= subject_length,
        "BLASTN HSP outside its sequences"
    );
    let context = if hsp.query_frame < 0 {
        reverse_complement(&query[query.len() - qe..query.len() - qs])
    } else {
        query[qs..qe].to_ascii_uppercase()
    };
    let subject = subject_letters(ss, se);
    let ungapped = [GapEditOp::Sub((qe - qs) as u32)];
    let ops = hsp.gap_info.as_deref().unwrap_or(&ungapped);
    // The rows live in a `PairwiseHit` while the reports are written: their exact size, not
    // what growing leaves (a third more on the alignments of the BLASTN hit records). A
    // reused buffer keeps its capacity.
    let columns = ops.iter().map(|op| op.num() as usize).sum();
    query_row.reserve_exact(columns);
    subject_row.reserve_exact(columns);
    // Offsets into the two segments.
    let (mut q, mut s) = (0, 0);
    for op in ops {
        let n = op.num() as usize;
        let (take_query, take_subject) = match op {
            GapEditOp::Sub(_) => (true, true),
            GapEditOp::Del(_) => (false, true),
            GapEditOp::Ins(_) => (true, false),
        };
        if take_query {
            query_row.extend_from_slice(
                context
                    .get(q..q + n)
                    .context("BLASTN edit script past the query")?,
            );
            q += n;
        } else {
            query_row.extend(repeat_n(b'-', n));
        }
        if take_subject {
            subject_row.extend_from_slice(
                subject
                    .get(s..s + n)
                    .context("BLASTN edit script past the subject")?,
            );
            s += n;
        } else {
            subject_row.extend(repeat_n(b'-', n));
        }
    }
    ensure!(
        q == context.len() && s == subject.len(),
        "BLASTN edit script does not span its HSP"
    );
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:878-892
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:825-842
// A mask that covers the whole query (its searched interval, `kTarget`) is not reported,
// so it is not shown. The masks of the search lie in the searched letters (`query_length`
// of them); `Map` moves them by the start of the query's interval (`offset`, 0 without
// `-query_loc`) into record coordinates, as the report's HSPs.
pub(crate) fn shown_query_masks(
    masks: &[MaskedInterval],
    query_length: usize,
    offset: usize,
) -> Vec<MaskedInterval> {
    merged(
        &masks
            .iter()
            .filter(|mask| !(mask.start == 0 && mask.end >= query_length))
            .map(|mask| MaskedInterval::new(mask.start + offset, mask.end + offset))
            .collect::<Vec<_>>(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1588-1600
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqvec.cpp:109-125
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqdb.cpp:121-133
// The subject masks shown with the HSPs of one subject (`GetMasks` of the subject set's
// `CSeqVecSeqInfoSrc` or of the database's `CSeqDbSeqInfoSrc`): each mask `[start, end)` that
// holds a letter (a stored database mask `TSeqRange(first, second-1)` of none is empty) and meets
// the closed range `[offset, end]` of an HSP in the searched letters (the end is one past the
// HSP), in the order of `masks`, merged for `lowercase_masked`, which marks each position once
// (the merged lowercase of a subject set's record is kept as it is). With `-subject_loc` the
// record's coordinates and the searched letters' differ by the interval's start, so a region of
// the record is kept by where it would lie in the interval, and shown where it lies in the
// record.
fn shown_subject_masks(masks: &[MaskedInterval], hsps: &[BlastnHsp]) -> Vec<MaskedInterval> {
    merged(&subject_mask_regions(masks, hsps).collect::<Vec<_>>())
}

/// The masks of `masks` that `GetSubjectMasks` gives with the HSPs `hsps` (`shown_subject_masks`
/// before the merge), in their order.
fn subject_mask_regions<'m>(
    masks: &'m [MaskedInterval],
    hsps: &'m [BlastnHsp],
) -> impl Iterator<Item = MaskedInterval> + 'm {
    masks
        .iter()
        .filter(|mask| {
            mask.start < mask.end
                && hsps.iter().any(|hsp| {
                    mask.start <= hsp.internal_s_end_0 && hsp.internal_s_offset_0 < mask.end
                })
        })
        .cloned()
}

/// The masked positions as disjoint intervals in order, for `lowercase_masked`.
fn merged(masks: &[MaskedInterval]) -> Vec<MaskedInterval> {
    let mut sorted = masks.to_vec();
    sorted.sort_by_key(|mask| mask.start);
    let mut merged: Vec<MaskedInterval> = Vec::with_capacity(sorted.len());
    for mask in sorted {
        match merged.last_mut() {
            Some(last) if mask.start <= last.end => last.end = last.end.max(mask.end),
            _ => merged.push(mask),
        }
    }
    merged
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2495-2521
// The row's residue positions start at `first` (1-based) and move by `step`; a residue
// whose position lies in a mask is shown in lowercase. `masks` are disjoint and in order
// (`merged`), so the masks that overlap the row's span are found by binary search.
fn lowercase_masked(row: &mut [u8], first: usize, step: isize, masks: &[MaskedInterval]) {
    let residues = row.iter().filter(|&&residue| residue != b'-').count();
    if masks.is_empty() || residues == 0 {
        return;
    }
    let last = first as isize + step * (residues as isize - 1);
    let (low, high) = (
        (first as isize).min(last) - 1,
        (first as isize).max(last) - 1,
    );
    let (low, high) = (low as usize, high as usize);
    let candidates = &masks[masks.partition_point(|mask| mask.end <= low)
        ..masks.partition_point(|mask| mask.start <= high)];
    if candidates.is_empty() {
        return;
    }
    let mut position = first as isize;
    for residue in row.iter_mut().filter(|residue| **residue != b'-') {
        let zero_based = (position - 1) as usize;
        if candidates
            .iter()
            .any(|mask| mask.start <= zero_based && zero_based < mask.end)
        {
            residue.make_ascii_lowercase();
        }
        position += step;
    }
}

/// The ungapped Karlin-Altschul block of a blastn query, from its composition; `None`
/// when NCBI cannot compute it and marks the query invalid.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2792
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2004-2013
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:346-350
/// The residues are counted in BLASTNA without `N` and the gap, normalized, and paired
/// with the standard composition (A, C, G, T at 0.25; `Blast_ResFreqStdComp`) through the
/// BLASTNA matrix (`BlastScoreFreqCalc`, blast_stat.c:2151-2205). The minus strand of a
/// query (its reverse complement) has the same residues in another order, so its block
/// can differ in the last bits; the search computes one per strand (`scoring.rs`). A
/// query of `N` only has no counted residue, so every score probability is 0 and the
/// calculation fails.
#[cfg(test)]
pub(crate) fn query_ungapped_karlin(
    query: &[u8],
    reward: i32,
    penalty: i32,
) -> Option<KarlinParams> {
    query_ungapped_karlin_blastna(
        &encode_iupac_to_blastna(&query.to_ascii_uppercase()),
        reward,
        penalty,
    )
}

/// `query_ungapped_karlin` of a context already in BLASTNA (the search holds its contexts
/// only in that encoding).
pub(crate) fn query_ungapped_karlin_blastna(
    blastna: &[u8],
    reward: i32,
    penalty: i32,
) -> Option<KarlinParams> {
    const BLASTNA_N: usize = 14;
    const BLASTNA_GAP: usize = 15;
    let mut counts = [0f64; 16];
    for &code in blastna {
        counts[(code & 0x0f) as usize] += 1.0;
    }
    counts[BLASTNA_N] = 0.0;
    counts[BLASTNA_GAP] = 0.0;
    let sum: f64 = counts.iter().sum();
    if sum == 0.0 {
        return None;
    }
    let matrix = build_blastna_matrix(reward, penalty);
    // At most 16 residues x 4 bases: a fixed array in place of a vector grown for each query.
    let mut probabilities = [(0i32, 0f64); 64];
    let mut num_probabilities = 0usize;
    for (residue, &count) in counts.iter().enumerate() {
        if count == 0.0 {
            continue;
        }
        for base in 0..4 {
            probabilities[num_probabilities] = (matrix[residue * 16 + base], count / sum * 0.25);
            num_probabilities += 1;
        }
    }
    let sfp = score_freq_profile_from_probabilities(
        penalty.min(reward),
        reward.max(penalty),
        &probabilities[..num_probabilities],
    );
    compute_karlin_params_ungapped(&sfp).ok()
}

#[cfg(test)]
#[path = "pairwise_tests.rs"]
mod tests;
