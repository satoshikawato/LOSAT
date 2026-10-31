//! `-lcase_masking` with a `-subject` that has lower-case letters: the subject's word scan over
//! NCBI's scan ranges, and its masked letters in the outfmt 0 rows. The case of a subject changes
//! nothing else, so a subject without lower-case letters keeps the search's loops, and only this
//! scan knows about ranges. Everything here is out
//! of line and not generic over the kernels (no NCBI counterpart: code placement only), so that the
//! loops of the default path keep their code. It reads the cells and the buffer in the order of
//! NCBI's scanners; the buffer's rounds change no order (`x_seed_fused_compressed.rs` header).

use std::collections::HashMap;

use super::super::x_seed_fused::x_test_hit;
use super::BlastpRuntimeLookup;
use crate::algorithm::tblastn::archive::lowercase_runs;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::X_BLASTP_SOFT;
use crate::algorithm::tblastx::lookup::{
    BlastAaLookupTable, XCell, XHit, PV_ARRAY_BTS, PV_ARRAY_MASK,
};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::pairwise::PairwiseHit;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:817-826
// ```c++
//             s_SeqLoc2MaskedSubjRanges(masks, &*range, length,  masked_ranges);
//             if ( !masked_ranges.empty() ) {
//                 ...
//             	subj->length = length;
//                 BlastSeqBlkSetSeqRanges(subj, (SSeqRange*) masked_ranges.get_data(),
//                                     static_cast<Uint4>(masked_ranges.size()) + 1, true, eSoftSubjMasking);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:211-213
// ```c
//     // Fill out the boundary of the sequence to compliment the masks
//     tmp[0].left = 0;
//     tmp[num_seq_ranges - 1].right = seq_blk->length;
// ```
/// The scan ranges of `letters` (a subject as it is searched, cut to `-subject_loc`): `(left,
/// right)`, one more than its lower-case runs. The array of pairs is read one place along, so each
/// range runs from the last letter of a run (inclusive) to the first letter of the next. Empty
/// without lower-case letters.
fn scan_ranges(letters: &[u8]) -> Vec<(i32, i32)> {
    let runs = lowercase_runs(letters);
    if runs.is_empty() {
        return Vec::new();
    }
    let mut ranges = Vec::with_capacity(runs.len() + 1);
    let mut left = 0;
    for (from, to) in runs {
        ranges.push((left, from));
        left = to;
    }
    ranges.push((left, letters.len() as i32));
    ranges
}

/// [`scan_ranges`] of each subject of a FASTA set.
#[cold]
#[inline(never)]
pub(super) fn x_subject_scan_ranges(records: &[FastaRecord]) -> Vec<Vec<(i32, i32)>> {
    records
        .iter()
        .map(|record| scan_ranges(record.seq()))
        .collect()
}

/// Makes `FastaRecord::cut` keep each record's whole letters (`archive_letters`) as well for a
/// `-subject_loc`: the outfmt 0 rows show the runs of the whole record ([`x_masked_subject_rows`]).
pub(super) fn x_keep_whole_letters(records: &mut [FastaRecord]) {
    for record in records {
        record.source.get_or_insert_with(Default::default);
    }
}

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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:75-127 (`s_BlastAaScanSubject`)
// ```c
//     while (s_DetermineScanningOffsets(subject, word_length, word_length, s_range)) {
//     s_first=subject->sequence + s_range[1];
//     s_last=subject->sequence + s_range[2];
//     /* prime the index */
//     index = ComputeTableIndex(word_length - 1, lookup->charsize, s_first);
//     for (s = s_first; s <= s_last; s++) {
//         index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
//         if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
// ```
/// The scan of the standard table over `ranges`, from `range` on: each hit in the order of NCBI's
/// buffer. A clamped first range can hold a word that runs past the sentinel, whose letters are read
/// as sentinels (as `x_seed_fused::x_scan_cells`).
fn scan_standard<K: XCell>(
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    ranges: &[(i32, i32)],
    range: &mut [i32; 3],
    visit: &mut dyn FnMut(OffsetPair),
) {
    let word_length = lookup.word_length;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let byte_at = |i: usize| subject.get(i).copied().unwrap_or(0) as usize;
    while range[1] <= range[2] {
        let mut index = 0usize;
        for i in 0..word_length as usize - 1 {
            index = (index << charsize) | byte_at(range[1] as usize + i);
        }
        for s in range[1]..=range[2] {
            let new_char = byte_at(s as usize + word_length as usize - 1);
            index = ((index << charsize) | new_char) & mask;
            if lookup.pv[index >> PV_ARRAY_BTS] & (1u32 << (index & PV_ARRAY_MASK)) == 0 {
                continue;
            }
            for hit in lookup.backbone[index].hits(&lookup.overflow) {
                visit(OffsetPair {
                    q_off: hit.q() as u32,
                    s_off: s as u32,
                });
            }
        }
        range[1] = range[2] + 1;
        while range[1] > range[2] {
            range[0] += 1;
            let Some(&(left, right)) = ranges.get(range[0] as usize) else {
                return;
            };
            range[1] = left;
            range[2] = right - word_length;
        }
    }
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
/// The word finders' scan over `ranges`: each word hit in NCBI's order. `clamps_scan_end` is the
/// two-hit search's raise of the first range's end, which scans one word of a first range too short
/// for it; the one-hit search scans nothing at all then, later ranges included.
#[cold]
#[inline(never)]
pub(super) fn x_for_each_hit(
    lookup: &BlastpRuntimeLookup,
    subject: &[u8],
    ranges: &[(i32, i32)],
    offset_pairs: &mut [OffsetPair],
    clamps_scan_end: bool,
    visit: &mut dyn FnMut(OffsetPair),
) {
    let word_size = lookup.word_length();
    let mut scan_range = [0, ranges[0].0, ranges[0].1 - word_size];
    if clamps_scan_end && scan_range[2] < scan_range[1] {
        scan_range[2] = scan_range[1];
    }
    let array_size =
        i32::try_from(offset_pairs.len()).expect("NCBI BLAST offset pair array size must fit Int4");
    match lookup {
        BlastpRuntimeLookup::Aa(table) => {
            scan_standard(table, subject, ranges, &mut scan_range, visit)
        }
        BlastpRuntimeLookup::Aa16(table) => {
            scan_standard(table, subject, ranges, &mut scan_range, visit)
        }
        BlastpRuntimeLookup::Compressed(table) => {
            while scan_range[1] <= scan_range[2] {
                let hits =
                    table.scan_subject(subject, ranges, offset_pairs, array_size, &mut scan_range);
                offset_pairs[..hits as usize]
                    .iter()
                    .for_each(|&pair| visit(pair));
            }
        }
        BlastpRuntimeLookup::CompressedFlat(table) => {
            while scan_range[1] <= scan_range[2] {
                let hits = table.scan_subject::<X_BLASTP_SOFT>(
                    subject,
                    ranges,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                offset_pairs[..hits as usize]
                    .iter()
                    .for_each(|&pair| visit(pair));
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-552 (quoted at
// `x_seed_fused::x_test_hit`)
/// The two-hit test of each hit of [`x_for_each_hit`], in scan order; `candidate(query_offset,
/// subject_offset, last_hit, diff, entry)` as for `x_seed_fused::x_scan_two_hit`.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_scan_two_hit_ranges<F>(
    lookup: &BlastpRuntimeLookup,
    subject: &[u8],
    ranges: &[(i32, i32)],
    offset_pairs: &mut [OffsetPair],
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    mut candidate: F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_for_each_hit(lookup, subject, ranges, offset_pairs, true, &mut |pair| {
        let query_offset = pair.q_off as i32;
        let subject_offset = pair.s_off as i32;
        let entry = &mut diag_array[((query_offset - subject_offset) & diag_mask) as usize];
        if let Some((last_hit, diff)) =
            x_test_hit(entry, subject_offset, diag_offset, window, word_size)
        {
            candidate(query_offset, subject_offset, last_hit, diff, entry);
        }
    });
}

// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1588-1600
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
// NCBI reference: c++/src/algo/blast/api/seqinfosrc_seqvec.cpp:115-122
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
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:2518-2521
// ```c++
//                     } else if (m_SeqLocChar==eLowerCase){
//                         actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
//                     }
// ```
/// The outfmt 0 copies of `hits` with each subject row's masked letters in lower case. A query's
/// hits with a subject show the lower-case runs of the subject's record (record coordinates) that
/// meet the closed range `[offset, end]` of one of their HSPs in the searched letters (the end is
/// one past the HSP): with `-subject_loc` the two coordinates differ by the interval's start, so a
/// run is kept by where it would lie in the interval and shown where it lies in the record. A `-db`
/// subject has no records here, and its rows stay as they are.
#[cold]
#[inline(never)]
pub(super) fn x_masked_subject_rows(
    hits: &[PairwiseHit],
    subject_records: &[FastaRecord],
    placements: &Placements,
) -> Vec<PairwiseHit> {
    let mut hsp_ranges: HashMap<(u32, u32), Vec<(i32, i32)>> = HashMap::new();
    for hit in hits {
        let shift = placements.offset(hit.hit.s_idx as usize);
        hsp_ranges
            .entry((hit.hit.q_idx, hit.hit.s_idx))
            .or_default()
            .push((
                (hit.hit.s_start - 1 - shift) as i32,
                (hit.hit.s_end - shift) as i32,
            ));
    }
    let mut hits = hits.to_vec();
    for hit in &mut hits {
        let Some(record) = subject_records.get(hit.hit.s_idx as usize) else {
            continue;
        };
        let Some(row) = hit.subject_seq.as_mut() else {
            continue;
        };
        let ranges = &hsp_ranges[&(hit.hit.q_idx, hit.hit.s_idx)];
        let runs: Vec<(i32, i32)> = lowercase_runs(record.archive_letters())
            .into_iter()
            .filter(|&(from, to)| {
                ranges
                    .iter()
                    .any(|&(offset, end)| from <= end && offset <= to)
            })
            .collect();
        if runs.is_empty() {
            continue;
        }
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
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:211-213 (quoted at
    // `scan_ranges`): the left end is the last letter of the run before, included.
    #[test]
    fn scan_ranges_run_from_the_last_lower_case_letter() {
        assert!(scan_ranges(b"MKV").is_empty());
        assert_eq!(scan_ranges(b"MKvvA"), [(0, 2), (3, 5)]);
        assert_eq!(scan_ranges(b"mKVvA"), [(0, 0), (0, 3), (3, 5)]);
        assert_eq!(scan_ranges(b"MKv"), [(0, 2), (2, 3)]);
        assert_eq!(scan_ranges(b"mkv"), [(0, 0), (2, 3)]);
    }
}
