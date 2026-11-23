//! `-lcase_masking` with a `-subject` that has lower-case letters: the subject's word scan over
//! NCBI's scan ranges, and its masked letters in the outfmt 0 rows. The case of a subject changes
//! nothing else, so a subject without lower-case letters keeps the search's loops, and only the
//! reference word finder knows about ranges (`seed::word_finder_over`). Everything here is out of
//! line (no NCBI counterpart: code placement only), so that the loops of the default path keep
//! their code. A `-db` subject has no case (the database stores the letters in upper case).

use anyhow::Result;

use super::args::ResolvedOptions;
use super::input::FastaRecord;
use super::parameters::ContextParameters;
use super::query_setup::PreparedQueryBatch;
use super::seed::{
    word_finder_over, x_compressed_scan, x_compressed_word_length, x_scan_bytes, BxLookup,
    Diagonals, InitHsp, Lookup,
};
use super::subjects::Subjects;
use crate::algorithm::tblastx::blast_aascan::s_blast_aa_scan_subject;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
use crate::report::pairwise::PairwiseHit;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:812-826
/// The scan ranges of subject `oid` when it is a record of `-subject` with lower-case runs
/// (`FastaRecord::lowercase_masks`, read with `-lcase_masking` only); `None` for any other subject,
/// whose scan is that of the whole subject.
#[inline]
pub(super) fn x_subject_scan_ranges(subjects: Subjects<'_>, oid: usize) -> Option<Vec<(i32, i32)>> {
    match subjects {
        Subjects::Fasta(records) if !records[oid].lowercase_masks.is_empty() => {
            Some(scan_ranges(&records[oid]))
        }
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:211-213
/// The scan ranges of `record`: `(left, right)`, one more than its lower-case runs. NCBI reads the
/// array of runs one place along, so each range runs from the last letter of a run (inclusive) to
/// the first letter of the next (as `blastp`'s `x_soft_subject::scan_ranges`).
#[cold]
#[inline(never)]
fn scan_ranges(record: &FastaRecord) -> Vec<(i32, i32)> {
    let mut ranges = Vec::with_capacity(record.lowercase_masks.len() + 1);
    let mut left = 0;
    for &(from, to) in &record.lowercase_masks {
        ranges.push((left, from));
        left = to;
    }
    ranges.push((left, record.sequence.len() as i32));
    ranges
}

/// The word finder of a subject with `ranges` (`seed::word_finder_over`), with each batch of word
/// hits added to `seeds` when it is given (the tests).
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn x_word_finder<S: super::seed::BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    ranges: &[(i32, i32)],
    mut seeds: Option<&mut Vec<(i32, i32)>>,
) -> Result<Vec<InitHsp>> {
    word_finder_over(
        scoring,
        batch,
        parameters,
        options,
        subject,
        |each| match lookup {
            BxLookup::Wide(lookup) => x_scan(lookup, subject, ranges, options.window_size, each),
            BxLookup::Narrow(lookup) => {
                lookup.x_scan(subject, 0, Some(ranges), options.window_size, each)
            }
        },
        diagonals,
        |pairs| {
            if let Some(seeds) = seeds.as_deref_mut() {
                seeds.extend_from_slice(pairs);
            }
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-505 (`s_BlastAaWordFinder_TwoHit`)
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-769 (`s_BlastAaWordFinder_OneHit`)
// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-59
/// The word finders' scan of `subject` (the residues and the sentinel) over `ranges`: `visit` with
/// each batch of the scanner's pairs, in NCBI's order. The two-hit search (`window` above 0) raises
/// the end of a first range too short for a word and scans that word; the one-hit search scans
/// nothing at all then, later ranges included. NCBI cuts a subject longer than `MAX_DBSEQ_LEN`
/// into chunks, each with ranges of its own (`blastp`'s `x_soft_subject::chunk_scans`); BLASTX
/// refuses such a subject (`search.rs`), so a subject here is one chunk.
#[cold]
#[inline(never)]
pub(super) fn x_scan(
    lookup: &Lookup,
    subject: &[u8],
    ranges: &[(i32, i32)],
    window: i32,
    visit: &mut dyn FnMut(&[(i32, i32)]),
) {
    let word = match lookup {
        Lookup::Standard(table, _) => table.word_length,
        Lookup::Compressed(table) => x_compressed_word_length(table),
    };
    let mut range = [0, ranges[0].0, ranges[0].1 - word];
    if window > 0 && range[2] < range[1] {
        range[2] = range[1];
    }
    let mut pad = Vec::new();
    let subject = x_scan_bytes(subject, range[2], word, &mut pad);
    let mut pairs = vec![OffsetPair::default(); lookup.capacity()];
    let capacity = pairs.len() as i32;
    let mut seeds = Vec::new();
    while range[1] <= range[2] {
        seeds.clear();
        let n = match lookup {
            Lookup::Standard(table, _) => {
                s_blast_aa_scan_subject(table, subject, ranges, &mut pairs, capacity, &mut range)
            }
            Lookup::Compressed(table) => {
                x_compressed_scan(table, subject, ranges, &mut pairs, capacity, &mut range)
            }
        };
        for pair in &pairs[..n as usize] {
            seeds.push((pair.q_off as i32, pair.s_off as i32));
        }
        visit(&seeds);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1588-1600
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqvec.cpp:115-122
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4302-4306 (`x_FillLocList`)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2518-2521
/// The subject rows' masks of `hits` (the outfmt 0 rows' `HspMasks`): NCBI keeps the subject
/// masks `kept` of the query and the subject (inclusive, in the subject's coordinates; the
/// query-anchored reports read the same lists): the lower-case runs of the subject's record that
/// meet one of the HSP ranges of the pair ([`x_subject_masks`], `-subject` read with
/// `-lcase_masking`), or the stored masks of a search with `-db_soft_mask` or `-db_hard_mask`
/// (`ReportedSubjects::shown_db_masks`). Each row places the kept masks that meet its own range
/// at its columns (`x_FillLocList`, `query_anchored::fill_row_loc_list`), which `x_OutputSeq`
/// puts in lower case. The rows come back so, as the hit records that the library API and LOSAT
/// Web hand out show them (`HspMasks::lower_row`); the lists are not kept
/// (`query_anchored::output_seq`): one `HspMasks` serves every row, except under `-html`, where
/// each subject row's go to its hit's entry of `html_masks` (by the place in `hits`,
/// `HspMasks::copy_row_from`). With `-subject_loc` the runs lie where they lie in the record, as
/// the hits' positions.
#[cold]
#[inline(never)]
pub(super) fn x_set_subject_row_masks(
    hits: &mut [PairwiseHit],
    kept: &std::collections::HashMap<(u32, u32), Vec<(i32, i32)>>,
    mut html_masks: Option<&mut Vec<crate::report::query_anchored::HspMasks>>,
) {
    let mut masks = crate::report::query_anchored::HspMasks::default();
    for (index, hit) in hits.iter_mut().enumerate() {
        let Some(runs) = kept.get(&(hit.hit.q_idx, hit.hit.s_idx)) else {
            continue;
        };
        // A row without letters (an ungapped HSP of length 0) has no mask.
        let Some(row) = hit.subject_seq.as_mut().filter(|row| !row.is_empty()) else {
            continue;
        };
        let mut letters = std::mem::take(row).into_bytes();
        // A stored `-` is a residue of the subject (`PairwiseHit::subject_stored_gaps`, in order).
        let stored = &hit.subject_stored_gaps;
        masks.set_row(
            1,
            0,
            &letters,
            |column| letters[column] != b'-' || stored.binary_search(&(column as u32)).is_ok(),
            (hit.hit.s_start as u32 - 1, hit.hit.s_end as u32 - 1),
            true,
            1,
            runs.iter().map(|&(from, to)| (from as u32, to as u32, 0)),
        );
        masks.lower_row(1, &mut letters);
        *row = String::from_utf8(letters).expect("ASCII alignment row");
        if let Some(entry) = html_masks
            .as_deref_mut()
            .and_then(|table| table.get_mut(index))
        {
            entry.copy_row_from(1, &masks);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1588-1600
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqvec.cpp:115-122
/// The subject masks of each query's hits (`GetSubjectMasks` of the query's results, by (query,
/// subject)): the lower-case runs of the subject's record (record coordinates, inclusive) that
/// meet the closed range `[offset, end]` of one of the query's HSPs with it in the searched letters
/// (the end one past the HSP): with `-subject_loc` the two coordinates differ by the interval's
/// start, so a run is kept by where it would lie in the interval and shown where it lies in the
/// record (BLASTP's `x_soft_subject::x_subject_masks`). Only records of `-subject` read with
/// `-lcase_masking` have runs. Outfmt 0's rows ([`x_set_subject_row_masks`]) and the
/// query-anchored reports (`report::query_anchored`) place them at their columns.
#[cold]
#[inline(never)]
pub(super) fn x_subject_masks(
    hits: &[PairwiseHit],
    subjects: &[FastaRecord],
) -> std::collections::HashMap<(u32, u32), Vec<(i32, i32)>> {
    let mut hsp_ranges: std::collections::HashMap<(u32, u32), Vec<(i32, i32)>> =
        std::collections::HashMap::new();
    for hit in hits.iter() {
        let Some(record) = subjects.get(hit.hit.s_idx as usize) else {
            continue;
        };
        let shift = record.placement().offset as i32;
        hsp_ranges
            .entry((hit.hit.q_idx, hit.hit.s_idx))
            .or_default()
            .push((
                hit.hit.s_start as i32 - 1 - shift,
                hit.hit.s_end as i32 - shift,
            ));
    }
    let mut record_runs: std::collections::HashMap<u32, Vec<(i32, i32)>> =
        std::collections::HashMap::new();
    hsp_ranges
        .into_iter()
        .filter_map(|(key, ranges)| {
            let record = subjects.get(key.1 as usize)?;
            let runs: Vec<(i32, i32)> = record_runs
                .entry(key.1)
                .or_insert_with(|| record.archive_lowercase_runs())
                .iter()
                .copied()
                .filter(|&(from, to)| {
                    ranges
                        .iter()
                        .any(|&(offset, end)| from <= end && offset <= to)
                })
                .collect();
            (!runs.is_empty()).then_some((key, runs))
        })
        .collect()
}

#[cfg(test)]
#[path = "x_soft_subject_tests.rs"]
mod tests;
