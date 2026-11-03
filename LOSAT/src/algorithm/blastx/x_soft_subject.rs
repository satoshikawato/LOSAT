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
    word_finder_over, x_compressed_scan, x_compressed_word_length, x_scan_bytes, Diagonals,
    InitHsp, Lookup,
};
use super::subjects::Subjects;
use crate::algorithm::tblastx::blast_aascan::s_blast_aa_scan_subject;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
use crate::report::pairwise::PairwiseHit;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:813-826
// ```c++
//         /* Set the lower case mask, if it exists */
//         if (subjects.GetMask(i).NotEmpty()) {
//             ...
//             s_SeqLoc2MaskedSubjRanges(masks, &*range, length,  masked_ranges);
//             if ( !masked_ranges.empty() ) {
//                 ...
//             	subj->length = length;
//                 BlastSeqBlkSetSeqRanges(subj, (SSeqRange*) masked_ranges.get_data(),
//                                     static_cast<Uint4>(masked_ranges.size()) + 1, true, eSoftSubjMasking);
// ```
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
// ```c
//     // Fill out the boundary of the sequence to compliment the masks
//     tmp[0].left = 0;
//     tmp[num_seq_ranges - 1].right = seq_blk->length;
// ```
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
/// hits added to `seeds` when the diagnostic stages trace them.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn x_word_finder(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    ranges: &[(i32, i32)],
    mut seeds: Option<&mut Vec<(i32, i32)>>,
) -> Result<Vec<InitHsp>> {
    word_finder_over(
        batch,
        parameters,
        options,
        subject,
        lookup,
        diagonals,
        Some(ranges),
        |pairs| {
            if let Some(seeds) = seeds.as_deref_mut() {
                seeds.extend_from_slice(pairs);
            }
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-505 (`s_BlastAaWordFinder_TwoHit`)
// ```c
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
//     scan_range[2] = subject->seq_ranges[0].right - wordsize;
//
//     if (scan_range[2] < scan_range[1])
//         scan_range[2] = scan_range[1];
//
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
//         hits = scansub(lookup_wrap, subject,
//                                   offset_pairs, array_size, scan_range);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-769 (`s_BlastAaWordFinder_OneHit`)
// ```c
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
//     scan_range[2] = subject->seq_ranges[0].right - wordsize;
//
//     while (scan_range[1] <= scan_range[2]) {
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-59
// ```c
//     while (range[1] > range[2]) {
//         range[0]++;
//         if (range[0] >= (Int4)subject->num_seq_ranges) {
//             return FALSE;
//         }
//         range[1] = subject->seq_ranges[range[0]].left + word_length - lut_word_length;
//         range[2] = subject->seq_ranges[range[0]].right - lut_word_length;
//     }
//     return TRUE;
// ```
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
// ```c++
//         // Union subject sequence ranges
//         vector <TSeqRange> ranges;
//         for (int i=0; i<hsp_list->hspcnt; i++) {
//             const BlastHSP* hsp = hsp_list->hsp_array[i];
//             TSeqRange rg;
//             rg.SetFrom(hsp->subject.offset);
//             rg.SetTo(hsp->subject.end);
//             ranges.push_back(rg);
//         }
//
//         // Extract subject masks
//         TMaskedSubjRegions masks;
//         if (!ranges.empty() && seqinfo_src->GetMasks(kOid, ranges, masks)) {
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqvec.cpp:115-122
// ```c++
//     for (size_t ir=0; ir< target_ranges.size(); ir++) {
//         if (target_ranges[ir] != TSeqRange::GetEmpty() &&
//            loc.IntersectingWith(target_ranges[ir])) {
//            CRef<CSeqLocInfo> sli(new CSeqLocInfo(interval, frame));
//            retval.push_back(sli);
//            return;
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4302-4304 (`x_FillLocList`)
// ```c++
//             if(interval.GetId().Match(m_AV->GetSeqId(i)) &&
//                m_AV->GetSeqRange(i).IntersectingWith(loc_range)){
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2518-2521
// ```c++
//                     } else if (m_SeqLocChar==eLowerCase){
//                         actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
//                     }
// ```
/// `hits` (the outfmt 0 rows) with each subject row's masked letters in lower case. NCBI keeps the
/// lower-case runs of the subject's record that meet one of the HSP ranges of the query and the
/// subject (`[offset, end]`, the end one past the HSP), then each row takes the kept runs that meet
/// its own range. A run that holds a letter of a row meets that row's HSP, so the rule is: every
/// letter of a row that lies in a run is in lower case (a run that only touches the end one past the
/// HSP holds no letter of the row). The runs (`FastaRecord::lowercase_masks`) are sorted and apart.
/// Only records of `-subject` read with `-lcase_masking` have runs; the caller calls this when one
/// has.
///
/// With `-subject_loc` the HSP ranges are in the searched letters and the runs in the whole record
/// (`FastaRecord::archive_lowercase_runs`): a run is kept by where it would lie in the interval and
/// shown where it lies in the record ([`x_shifted_subject_rows`], BLASTP's
/// `x_soft_subject::x_masked_subject_rows`). An interval that starts at the record's start moves
/// nothing, and the rule above holds for its clipped runs.
#[cold]
#[inline(never)]
pub(super) fn x_masked_subject_rows(hits: &mut [PairwiseHit], subjects: &[FastaRecord]) {
    if subjects
        .iter()
        .any(|subject| subject.placement().offset > 0)
    {
        return x_shifted_subject_rows(hits, subjects);
    }
    for hit in hits.iter_mut() {
        let Some(record) = subjects.get(hit.hit.s_idx as usize) else {
            continue;
        };
        let runs = &record.lowercase_masks;
        let mut position = hit.hit.s_start as i32 - 1;
        // The first run that ends at the row's first letter or after it.
        let mut next = runs.partition_point(|&(_, to)| to < position);
        if runs
            .get(next)
            .is_none_or(|&(from, _)| from > hit.hit.s_end as i32 - 1)
        {
            continue;
        }
        let Some(row) = hit.subject_seq.as_mut() else {
            continue;
        };
        let mut bytes = std::mem::take(row).into_bytes();
        for (column, letter) in bytes.iter_mut().enumerate() {
            // A stored `-` is a residue of the subject (`PairwiseHit::subject_stored_gaps`).
            if *letter == b'-' && !hit.subject_stored_gaps.contains(&(column as u32)) {
                continue;
            }
            while runs.get(next).is_some_and(|&(_, to)| to < position) {
                next += 1;
            }
            if runs.get(next).is_some_and(|&(from, _)| from <= position) {
                letter.make_ascii_lowercase();
            }
            position += 1;
        }
        *row = String::from_utf8(bytes).expect("ASCII alignment row");
    }
}

/// [`x_masked_subject_rows`] of subjects cut to a `-subject_loc` interval that starts past the
/// record's start. A query's hits with a subject keep the runs of the whole record (record
/// coordinates) that meet the closed range `[offset, end]` of one of their HSPs in the searched
/// letters (the end one past the HSP); each row shows the kept runs where they lie in the record.
fn x_shifted_subject_rows(hits: &mut [PairwiseHit], subjects: &[FastaRecord]) {
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
    for hit in hits.iter_mut() {
        let Some(record) = subjects.get(hit.hit.s_idx as usize) else {
            continue;
        };
        let ranges = &hsp_ranges[&(hit.hit.q_idx, hit.hit.s_idx)];
        let runs: Vec<(i32, i32)> = record_runs
            .entry(hit.hit.s_idx)
            .or_insert_with(|| record.archive_lowercase_runs())
            .iter()
            .copied()
            .filter(|&(from, to)| {
                ranges
                    .iter()
                    .any(|&(offset, end)| from <= end && offset <= to)
            })
            .collect();
        if runs.is_empty() {
            continue;
        }
        let Some(row) = hit.subject_seq.as_mut() else {
            continue;
        };
        let mut bytes = std::mem::take(row).into_bytes();
        let mut position = hit.hit.s_start as i32 - 1;
        for (column, letter) in bytes.iter_mut().enumerate() {
            // A stored `-` is a residue of the subject (`PairwiseHit::subject_stored_gaps`).
            if *letter == b'-' && !hit.subject_stored_gaps.contains(&(column as u32)) {
                continue;
            }
            if runs
                .iter()
                .any(|&(from, to)| from <= position && position <= to)
            {
                letter.make_ascii_lowercase();
            }
            position += 1;
        }
        *row = String::from_utf8(bytes).expect("ASCII alignment row");
    }
}

#[cfg(test)]
mod tests {
    use super::super::input::parse_fasta;
    use super::super::{parameters, query_setup::prepare_queries};
    use super::*;

    /// A blastx search's lookup table, contexts and options for one nucleotide query of `letters`
    /// with `extra` arguments (the task), as `seed.rs`'s tests build them.
    fn search(
        letters: &[u8],
        extra: &[&str],
    ) -> (
        PreparedQueryBatch,
        Vec<ContextParameters>,
        ResolvedOptions,
        Lookup,
    ) {
        let fasta = [b">q\n".as_slice(), letters, b"\n"].concat();
        let q = parse_fasta(&fasta, false, false).unwrap();
        let mut args = vec![
            "losat", "blastx", "-query", "query", "-subject", "subject", "-seg", "no",
        ];
        args.extend_from_slice(extra);
        let crate::cli::Commands::Blastx(args) =
            crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
                .unwrap()
                .command
        else {
            unreachable!()
        };
        let options = args.resolve().unwrap();
        let batch = prepare_queries(&q, &options).unwrap();
        let mut params = parameters::score_block(&batch);
        parameters::effective_lengths(&batch, &mut params, &options, 1_000_000, 1, None).unwrap();
        let lookup = Lookup::new(&batch, &params, &options).unwrap();
        (batch, params, options, lookup)
    }

    /// `length` letters of a fixed pseudo-random protein, encoded, with the sentinel after them.
    fn subject(length: usize) -> Vec<u8> {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut letters: Vec<u8> = (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                crate::utils::matrix::aa_char_to_ncbistdaa(
                    b"ARNDCQEGHILKMFPSTWYV"[(state % 20) as usize],
                )
            })
            .collect();
        letters.push(0);
        letters
    }

    /// The word hits of the word finders' scan of `subject` over `ranges` (the whole subject
    /// without).
    fn seeds(
        lookup: &Lookup,
        subject: &[u8],
        ranges: Option<&[(i32, i32)]>,
        window: i32,
    ) -> Vec<(i32, i32)> {
        let mut seeds = Vec::new();
        match ranges {
            None => lookup.scan(subject, subject.len() - 1, window, |pairs| {
                seeds.extend_from_slice(pairs)
            }),
            Some(ranges) => x_scan(lookup, subject, ranges, window, &mut |pairs| {
                seeds.extend_from_slice(pairs)
            }),
        }
        seeds
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-59 and
    // aa_ungapped.c:495-505, :762-769 (quoted at `x_scan`): a range `(left, right)` holds the words
    // that start at `left..=right - word`; a first range too short for a word holds its first word
    // in the two-hit search and stops the one-hit search's scan.
    #[test]
    fn the_scan_over_ranges_takes_the_words_of_the_ranges() {
        let query = b"ATGGCTCGTAACGATTGCCAAGAAGGACATATCCTGAAAATGTTCCCAAGCACCTGGTACGTCGCTCGTAACGATTGCCAAGAAGGACATATCCTGAAAATGTTCCCAAGCACCTGGTACGTCGCTCGTAACGATTGCCAAGAAGGACATATCCTG";
        for extra in [&[][..], &["-task", "blastx-fast"][..]] {
            let (_, _, options, lookup) = search(query, extra);
            let word = options.word_size;
            for length in [400usize, 2] {
                let subject = subject(length);
                let length = length as i32;
                for window in [0, 40] {
                    let whole = seeds(&lookup, &subject, None, window);
                    assert_eq!(
                        seeds(&lookup, &subject, Some(&[(0, length)]), window),
                        whole
                    );
                    let runs: &[&[(i32, i32)]] = if length == 2 {
                        &[&[(1, 1)], &[(0, 1)]]
                    } else {
                        &[
                            &[(10, 20)],
                            &[(0, 0), (50, 51)],
                            &[(2, 30), (100, 200)],
                            &[(length - 5, length - 1)],
                            &[(0, length - 1)],
                        ]
                    };
                    for &runs in runs {
                        let ranges = scan_ranges(&record(&vec![b'A'; length as usize], runs));
                        let short = ranges[0].1 - word < ranges[0].0;
                        let expected: Vec<(i32, i32)> = whole
                            .iter()
                            .copied()
                            .filter(|&(_, s): &(i32, i32)| {
                                (short && window > 0 && s == 0)
                                    || ranges
                                        .iter()
                                        .any(|&(left, right)| left <= s && s <= right - word)
                            })
                            .filter(|_| !(short && window == 0))
                            .collect();
                        assert_eq!(
                            seeds(&lookup, &subject, Some(&ranges), window),
                            expected,
                            "{extra:?} length {length} window {window} runs {runs:?}"
                        );
                    }
                }
            }
        }
    }

    // The reference word finder over one range that covers each subject is the word finder of the
    // whole subjects: the same hits, here and on the subjects after it (the same diagonal table).
    #[test]
    fn one_range_over_the_subject_is_the_whole_scan() {
        let query = b"ATGGCTCGTAACGATTGCCAAGAAGGACATATCCTGAAAATGTTCCCAAGCACCTGGTACGTCGCTCGTAACGATTGCCAAGAAGGACATATCCTGAAAATGTTCCCAAGCACCTGGTACGTCGCTCGTAACGATTGCCAAGAAGGACATATCCTG";
        for window in ["0", "40"] {
            let (batch, params, options, lookup) = search(query, &["-window_size", window]);
            let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
            let mut whole = Diagonals::new(length, options.window_size).unwrap();
            let mut ranged = Diagonals::new(length, options.window_size).unwrap();
            for letters in [300, 7, 300] {
                let subject = subject(letters);
                let expected = super::super::seed::word_finder(
                    &batch,
                    &params,
                    &options,
                    &subject,
                    &lookup,
                    &mut whole,
                    |_| {},
                )
                .unwrap();
                let mut seeds = Vec::new();
                let hits = x_word_finder(
                    &batch,
                    &params,
                    &options,
                    &subject,
                    &lookup,
                    &mut ranged,
                    &[(0, letters as i32)],
                    Some(&mut seeds),
                )
                .unwrap();
                assert!(!seeds.is_empty());
                assert_eq!(hits, expected, "window {window} letters {letters}");
            }
        }
    }

    fn record(letters: &[u8], runs: &[(i32, i32)]) -> FastaRecord {
        FastaRecord {
            internal_id: String::new(),
            title: String::new(),
            sequence: letters.to_vec(),
            lowercase_masks: runs.to_vec(),
            warnings: Vec::new(),
            source: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:211-213 (quoted at
    // `scan_ranges`): the left end is the last letter of the run before, included.
    #[test]
    fn scan_ranges_run_from_the_last_lower_case_letter() {
        assert_eq!(scan_ranges(&record(b"MKVVA", &[(2, 3)])), [(0, 2), (3, 5)]);
        assert_eq!(
            scan_ranges(&record(b"MKVVA", &[(0, 0), (3, 3)])),
            [(0, 0), (0, 3), (3, 5)]
        );
        assert_eq!(scan_ranges(&record(b"MKV", &[(2, 2)])), [(0, 2), (2, 3)]);
        assert_eq!(scan_ranges(&record(b"MKV", &[(0, 2)])), [(0, 0), (2, 3)]);
    }

    #[test]
    fn only_fasta_records_with_runs_have_ranges() {
        let records = [record(b"MKVVA", &[]), record(b"MKVVA", &[(1, 2)])];
        let subjects = Subjects::Fasta(&records);
        assert_eq!(x_subject_scan_ranges(subjects, 0), None);
        assert_eq!(
            x_subject_scan_ranges(subjects, 1),
            Some(vec![(0, 1), (2, 5)])
        );
        assert_eq!(x_subject_scan_ranges(Subjects::Unreadable, 0), None);
    }
}
