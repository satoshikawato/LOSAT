//! The value rules of an `Hsp` (`s_SeqAlignSetToXMLHsps`, blastxml_format.cpp:107-339): the
//! six-digit scores, the masked query row and the recount of identities and positives.

use super::matrix::FormattingMatrix;
use super::records::QueryMask;

/// `sprintf("%.*g", 6, x)` then `atof` (blastxml_format.cpp:131-136):
/// ```c++
///         //Print 6 significant digits for double values
///         char tmp[512];
///         sprintf(tmp,"%.*g", 6, bit_score );
///         bit_score = atof(tmp);
///         sprintf(tmp,"%.*g", 6, evalue );
///         evalue = atof(tmp);
/// ```
/// Both conversions are correctly rounded (glibc's `printf` and `strtod`; Rust's `{:e}` with a
/// precision and `parse`), ties to even on the exact binary value.
pub fn round6(value: f64) -> f64 {
    if !value.is_finite() || value == 0.0 {
        return value;
    }
    format!("{value:.5e}").parse().unwrap_or(value)
}

/// The query row of an HSP with its masked regions replaced by the mask letter
/// (`s_MaskQuerySeq`, blastfmtutil.cpp:446-521), literally.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:446-521
/// ```c++
///     vector<CRange<int> > masks_v;
///     int aln_stop = static_cast<int>(query_seq.size()) - 1;
///     ITERATE(ncbi::TMaskedQueryRegions, mask_iter, mask_info) {
///         if ((*mask_iter)->GetFrame() != query_frame)
///             continue;
///         int start =
///             alnvec.GetAlnPosFromSeqPos(0,
///                                        (*mask_iter)->GetInterval().GetFrom());
///         int stop =
///             alnvec.GetAlnPosFromSeqPos(0,
///                                        (*mask_iter)->GetInterval().GetTo());
///         // For negative frames, start and stop must be swapped.
///         ...
///         if (start >= 0) {
///             if (stop < 0)
///                 stop = aln_stop;
///             CRange<int>  range(start, stop);
///             masks_v.push_back(range);
///         }
///     }
///
///     sort(masks_v.begin(), masks_v.end(), SRangeStartSort());
///
///     // Mask the sequence
///     int mask_index = 0;
///     for (int seg_index = 0;
///          seg_index < (int) segs_v.size() && mask_index < (int) masks_v.size();
///          ++seg_index) {
///         ...
///         while (mask_index < (int) masks_v.size() &&
///                (mask_pos = max(seg_start, masks_v[mask_index].GetFrom()))
///                <= seg_stop) {
///             int mask_stop = min(seg_stop, masks_v[mask_index].GetTo());
///             // Mask the respective part of the sequence
///             for ( ; mask_pos <= mask_stop; ++mask_pos) {
///                 if(  query_seq[mask_pos] == '-' ) continue; // preserve gap
///                 ...
///             }
///             // Advance to the next mask if this mask is done with. Otherwise
///             // break out of the loop.
///             if (mask_pos < seg_stop)
///                 ++mask_index;
///             else
///                 break;
///         }
///     }
/// ```
/// `GetAlnPosFromSeqPos(0, position)` is -1 for a position outside the HSP's query range, so a
/// mask that starts before the HSP is dropped whole, and a mask that ends after it is cut at
/// the last column.
///
/// `row` is the unmasked query row, `segments` the Dense-seg segments as inclusive column
/// ranges (`GetAlnStart` / `GetAlnStop`), `first` the 0-based record position of the row's
/// first residue (query residues rise by one per residue: the protein and plus-strand rows),
/// `masks` the query's masked regions in record coordinates. The masked row is written to
/// `out`; `columns` is scratch.
#[allow(clippy::too_many_arguments)]
pub fn mask_query_row(
    row: &[u8],
    segments: &[(usize, usize)],
    first: i32,
    masks: &[QueryMask],
    query_frame: i32,
    mask_char: u8,
    out: &mut Vec<u8>,
    columns: &mut Vec<usize>,
) {
    out.clear();
    out.extend_from_slice(row);
    if !masks.iter().any(|mask| mask.frame == query_frame) {
        return;
    }
    // The column of each query residue of the row, in row order.
    columns.clear();
    columns.extend(
        row.iter()
            .enumerate()
            .filter(|&(_, &letter)| letter != b'-')
            .map(|(column, _)| column),
    );
    let aln_pos = |position: i32| -> i32 {
        let index = i64::from(position) - i64::from(first);
        if index < 0 || index >= columns.len() as i64 {
            -1
        } else {
            columns[index as usize] as i32
        }
    };
    let aln_stop = row.len() as i32 - 1;
    let mut ranges: Vec<(i32, i32)> = Vec::new();
    for mask in masks.iter().filter(|mask| mask.frame == query_frame) {
        let mut start = aln_pos(mask.from);
        let mut stop = aln_pos(mask.to);
        if query_frame < 0 {
            std::mem::swap(&mut start, &mut stop);
        }
        if start >= 0 {
            if stop < 0 {
                stop = aln_stop;
            }
            ranges.push((start, stop));
        }
    }
    ranges.sort_by_key(|&(start, _)| start);
    let mut mask_index = 0;
    for &(seg_start, seg_stop) in segments {
        if mask_index >= ranges.len() {
            break;
        }
        let (seg_start, seg_stop) = (seg_start as i32, seg_stop as i32);
        if seg_start > seg_stop {
            continue;
        }
        while mask_index < ranges.len() {
            let mut mask_pos = seg_start.max(ranges[mask_index].0);
            if mask_pos > seg_stop {
                break;
            }
            let mask_stop = seg_stop.min(ranges[mask_index].1);
            while mask_pos <= mask_stop {
                let letter = &mut out[mask_pos as usize];
                if *letter != b'-' {
                    *letter = mask_char;
                }
                mask_pos += 1;
            }
            if mask_pos < seg_stop {
                mask_index += 1;
            } else {
                break;
            }
        }
    }
}

/// The identities, positives and midline of an HSP, recounted from the unmasked rows
/// (blastxml_format.cpp:296-331):
/// ```c++
///         for (unsigned int i = 0; i < kMaxOffset; ++i) {
///             if (query_seq[i] == subject_seq[i]) {
///                 ++num_ident;
///                 ++num_positives;
///                 if (kIsBlastn)
///                     middle_seq[i] = '|';
///             } else if (matrix &&
///                        (*matrix)(query_seq[i], subject_seq[i]) > 0 &&
///                        !kIsBlastn) {
///                 ++num_positives;
///                 middle_seq[i] = kIsBlastn ? ' ' : '+';
///             } else {
///                 middle_seq[i] = ' ';
///             }
///         }
/// ```
/// `matrix` is `None` for BLASTN. The midline starts as a copy of the query row (an identity
/// keeps its letter). Returns `(identity, positive)`.
pub fn recount(
    qrow: &[u8],
    srow: &[u8],
    matrix: Option<&FormattingMatrix>,
    midline: &mut Vec<u8>,
) -> (i32, i32) {
    midline.clear();
    midline.extend_from_slice(qrow);
    let (mut identity, mut positive) = (0, 0);
    for (i, (&q, &s)) in qrow.iter().zip(srow).enumerate() {
        if q == s {
            identity += 1;
            positive += 1;
            if matrix.is_none() {
                midline[i] = b'|';
            }
        } else if matrix.is_some_and(|matrix| matrix.positive(q, s)) {
            positive += 1;
            midline[i] = b'+';
        } else {
            midline[i] = b' ';
        }
    }
    (identity, positive)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mask(from: i32, to: i32) -> QueryMask {
        QueryMask { from, to, frame: 0 }
    }

    fn masked(row: &[u8], segments: &[(usize, usize)], first: i32, masks: &[QueryMask]) -> String {
        let (mut out, mut columns) = (Vec::new(), Vec::new());
        mask_query_row(row, segments, first, masks, 0, b'X', &mut out, &mut columns);
        String::from_utf8(out).expect("ASCII")
    }

    #[test]
    fn round6_rounds_to_six_digits_and_ties_to_even() {
        assert_eq!(round6(255.7579999), 255.758);
        assert_eq!(round6(5.303258e-94), 5.30326e-94);
        assert_eq!(round6(448.0104), 448.01);
        // exact ties go to the even digit, as glibc's printf does
        assert_eq!(round6(1234565.0), 1234560.0);
        assert_eq!(round6(1234575.0), 1234580.0);
        assert_eq!(round6(0.0), 0.0);
        assert_eq!(round6(1e-5), 1e-5);
        assert_eq!(round6(10.0), 10.0);
        assert_eq!(round6(123456789.0), 123457000.0);
    }

    #[test]
    fn recount_counts_letters_and_positives() {
        let matrix = FormattingMatrix::new("BLOSUM62").unwrap();
        let mut midline = Vec::new();
        // identity, I/V positive, A/W not, J/V positive (XML only), gap columns, identity
        let (identity, positive) = recount(b"MIAJ-UK", b"MVWVK-K", Some(&matrix), &mut midline);
        assert_eq!(String::from_utf8(midline).unwrap(), "M+ +  K");
        assert_eq!((identity, positive), (2, 4));
    }

    #[test]
    fn blastn_recount_has_equal_identity_and_positive() {
        let mut midline = Vec::new();
        let (identity, positive) = recount(b"ACGT-AC", b"ACCTTAC", None, &mut midline);
        assert_eq!(String::from_utf8(midline).unwrap(), "|| | ||");
        assert_eq!((identity, positive), (5, 5));
    }

    #[test]
    fn a_mask_inside_the_hsp_replaces_the_letters() {
        // HSP of query positions 10..=17 (first = 10) with a gap column: columns 0-2, gap, 4-8.
        let row = b"ACD-EFGHI";
        let segments = [(0, 2), (3, 3), (4, 8)];
        // mask positions 11..=12 -> columns 1, 2
        assert_eq!(masked(row, &segments, 10, &[mask(11, 12)]), "AXX-EFGHI");
        // mask positions 12..=13 spans the gap column (kept as '-')
        assert_eq!(masked(row, &segments, 10, &[mask(12, 13)]), "ACX-XFGHI");
    }

    #[test]
    fn a_mask_starting_before_the_hsp_is_dropped_whole() {
        let row = b"ACDEFGHI";
        let segments = [(0, 7)];
        // the mask 8..=12 covers the HSP's first positions 10-12 but starts before it
        assert_eq!(masked(row, &segments, 10, &[mask(8, 12)]), "ACDEFGHI");
        // a mask starting inside and running past the end is cut at the last column
        assert_eq!(masked(row, &segments, 10, &[mask(15, 40)]), "ACDEFXXX");
        // a mask wholly outside is ignored
        assert_eq!(masked(row, &segments, 10, &[mask(30, 40)]), "ACDEFGHI");
    }

    #[test]
    fn masks_of_another_frame_are_ignored() {
        let (mut out, mut columns) = (Vec::new(), Vec::new());
        let masks = [QueryMask {
            from: 11,
            to: 12,
            frame: 1,
        }];
        mask_query_row(
            b"ACDEF",
            &[(0, 4)],
            10,
            &masks,
            0,
            b'X',
            &mut out,
            &mut columns,
        );
        assert_eq!(out, b"ACDEF");
        mask_query_row(
            b"ACDEF",
            &[(0, 4)],
            10,
            &masks,
            1,
            b'N',
            &mut out,
            &mut columns,
        );
        assert_eq!(out, b"ANNEF");
    }

    #[test]
    fn the_segment_loop_is_ncbis_literally() {
        // A mask that ends on a segment's second-last column does not advance (mask_pos ==
        // seg_stop), so the next segment starts the same mask again from its own first column:
        // max(seg_start, from) is then past mask_stop and nothing more is masked.
        let row = b"ABCDEFGH";
        let segments = [(0, 3), (4, 7)];
        assert_eq!(masked(row, &segments, 0, &[mask(1, 2)]), "AXXDEFGH");
        // Two masks in one segment: the first ends before the segment's last column.
        assert_eq!(
            masked(row, &segments, 0, &[mask(0, 0), mask(2, 3)]),
            "XBXXEFGH"
        );
        // A one-column segment: mask_pos = seg_stop + 1 after the masking, so `break`.
        let segments = [(0, 3), (4, 4), (5, 7)];
        assert_eq!(masked(row, &segments, 0, &[mask(4, 6)]), "ABCDXXXH");
        // Masks of the same start keep their relative order (a stable sort).
        assert_eq!(
            masked(row, &[(0, 7)], 0, &[mask(1, 2), mask(1, 4)]),
            "AXXXXFGH"
        );
    }
}
