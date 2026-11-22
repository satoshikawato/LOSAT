//! The value rules of an `Hsp` (`s_SeqAlignSetToXMLHsps`, blastxml_format.cpp:107-339): the
//! six-digit scores, the masked query row and the recount of identities and positives.

use super::matrix::FormattingMatrix;
use super::records::QueryMask;

/// `sprintf` with the format "%.*g" and the precision 6, then `atof` (blastxml_format.cpp:131-136):
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:119-124
/// Both conversions are correctly rounded (glibc's `printf` and `strtod`; Rust's `{:e}` with a
/// precision and `parse`), ties to even on the exact binary value.
pub fn round6(value: f64) -> f64 {
    if !value.is_finite() || value == 0.0 {
        return value;
    }
    format!("{value:.5e}").parse().unwrap_or(value)
}

/// What a masked query letter becomes: outfmt 5 writes one fixed letter (`X`, `N`;
/// `eSubstituteMask`... `s_MaskQuerySeq`'s `eXMask` / `eNMask` forms), XML2 the lower-case
/// letter (`eLowerCase`, blastxml2_format.cpp:279); `tolower` leaves `*` as it is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaskLetter {
    Fixed(u8),
    Lower,
}

/// The query row of an HSP with its masked regions replaced by the mask letter
/// (`s_MaskQuerySeq`, blastfmtutil.cpp:446-521), literally.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:446-521
/// `GetAlnPosFromSeqPos` of row 0 and `position` is -1 for a position outside the HSP's query
/// range, so a mask that starts before the HSP is dropped whole, and a mask that ends after it is
/// cut at the last column.
///
/// `row` is the unmasked query row, `segments` the Dense-seg segments as inclusive column
/// ranges (`GetAlnStart` / `GetAlnStop`), `masks` the query's masked regions in record
/// coordinates. `first` and `step` place the row's residues in the record: residue `j` (the
/// `j`-th letter of the row that is not a gap) lies at `first + j` for `step` 1 (a protein or
/// plus-strand nucleotide row); at the nucleotides `first + 3j ..= first + 3j + 2` for `step`
/// 3 (a translated query, plus frame); at `first - 3j - 2 ..= first - 3j` for `step` -3 (minus
/// frame, the row runs from the high end: `first` is the highest nucleotide). A position
/// maps to the column of the residue that holds it
/// (the delta is `seq_pos` less `start`, divided by `GetWidth` of the row, alnmap.cpp:593-595), and
/// to -1 outside the row's residues. The masked row is written to `out`; `columns` is scratch.
#[allow(clippy::too_many_arguments)]
pub fn mask_query_row(
    row: &[u8],
    segments: &[(usize, usize)],
    first: i32,
    step: i32,
    masks: &[QueryMask],
    query_frame: i32,
    mask_letter: MaskLetter,
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
    let width = i64::from(step.abs().max(1));
    let aln_pos = |position: i32| -> i32 {
        let distance = if step < 0 {
            i64::from(first) - i64::from(position)
        } else {
            i64::from(position) - i64::from(first)
        };
        if distance < 0 || distance / width >= columns.len() as i64 {
            -1
        } else {
            columns[(distance / width) as usize] as i32
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
    // NCBI's `std::sort` with `SRangeStartSort` is not stable, so the order of two masks of one
    // start is unspecified there; this one keeps their order. Harmless for BLASTP (the SEG
    // masks of a query are merged first, so no two start together).
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
                    *letter = match mask_letter {
                        MaskLetter::Fixed(letter) => letter,
                        MaskLetter::Lower => letter.to_ascii_lowercase(),
                    };
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:302-314
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

/// `s_GetTranslationFrame` (blastxml_format.cpp:85-97) on the 1-based endpoints of a translated
/// HSP's interval in a record (`start` the lower, `end` the higher) of `length` letters: the
/// reading frame by the record's own coordinates, whatever interval was searched.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:88-96
pub fn translation_frame(plus_strand: bool, start: i32, end: i32, length: i32) -> i32 {
    if plus_strand {
        (start - 1) % 3 + 1
    } else {
        -((length - end) % 3 + 1)
    }
}

#[cfg(test)]
#[path = "hsp_tests.rs"]
mod tests;
