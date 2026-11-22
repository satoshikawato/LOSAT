use super::*;

fn mask(from: i32, to: i32) -> QueryMask {
    QueryMask { from, to, frame: 0 }
}

fn masked(row: &[u8], segments: &[(usize, usize)], first: i32, masks: &[QueryMask]) -> String {
    let (mut out, mut columns) = (Vec::new(), Vec::new());
    mask_query_row(
        row,
        segments,
        first,
        1,
        masks,
        0,
        MaskLetter::Fixed(b'X'),
        &mut out,
        &mut columns,
    );
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
fn the_lower_case_mask_letter_keeps_a_star() {
    // XML2: `tolower`, which leaves `*` and the gap as they are (tblastx_db_gc HSP 14)
    let (mut out, mut columns) = (Vec::new(), Vec::new());
    mask_query_row(
        b"AC*-EF",
        &[(0, 5)],
        0,
        1,
        &[mask(0, 3)],
        0,
        MaskLetter::Lower,
        &mut out,
        &mut columns,
    );
    assert_eq!(out, b"ac*-eF");
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
        1,
        &masks,
        0,
        MaskLetter::Fixed(b'X'),
        &mut out,
        &mut columns,
    );
    assert_eq!(out, b"ACDEF");
    mask_query_row(
        b"ACDEF",
        &[(0, 4)],
        10,
        1,
        &masks,
        1,
        MaskLetter::Fixed(b'N'),
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
    // Masks that are given out of order are sorted by their start (`SRangeStartSort`; the
    // order of equal starts is unspecified in NCBI and not tested here).
    assert_eq!(
        masked(row, &[(0, 7)], 0, &[mask(5, 6), mask(1, 2)]),
        "AXXDEXXH"
    );
}

/// A translated query row: nucleotide masks reach the columns of the residues that hold
/// them (`GetAlnPosFromSeqPos`: `(seq_pos - start) / width`), the minus frame's row runs
/// from the highest nucleotide, and a mask whose first column in alignment order lies
/// before the HSP is dropped.
#[test]
fn nucleotide_masks_of_a_translated_row_reach_the_columns_of_their_codons() {
    let frame = |frame: i32, from: i32, to: i32| QueryMask { from, to, frame };
    let run = |row: &[u8], first: i32, step: i32, query_frame: i32, masks: &[QueryMask]| {
        let (mut out, mut columns) = (Vec::new(), Vec::new());
        let segments = [(0, row.len() - 1)];
        mask_query_row(
            row,
            &segments,
            first,
            step,
            masks,
            query_frame,
            MaskLetter::Fixed(b'X'),
            &mut out,
            &mut columns,
        );
        String::from_utf8(out).unwrap()
    };
    // Plus frame: residue j holds the nucleotides 30 + 3j ..= 32 + 3j.
    assert_eq!(run(b"ABCDE", 30, 3, 1, &[frame(1, 36, 41)]), "ABXXE");
    // a mask inside a codon masks that codon's residue
    assert_eq!(run(b"ABCDE", 30, 3, 1, &[frame(1, 37, 37)]), "ABXDE");
    // starting before the HSP: dropped whole; running past its end: cut at the last column
    assert_eq!(run(b"ABCDE", 30, 3, 1, &[frame(1, 28, 36)]), "ABCDE");
    assert_eq!(run(b"ABCDE", 30, 3, 1, &[frame(1, 40, 90)]), "ABCXX");
    // another frame's mask is not applied
    assert_eq!(run(b"ABCDE", 30, 3, 1, &[frame(2, 36, 41)]), "ABCDE");
    // Minus frame: the row runs from the highest nucleotide 47 down; residue j holds
    // 45 - 3j ..= 47 - 3j.
    assert_eq!(run(b"ABCDE", 47, -3, -1, &[frame(-1, 36, 41)]), "ABXXE");
    assert_eq!(run(b"ABCDE", 47, -3, -1, &[frame(-1, 30, 41)]), "ABXXX");
    // reaching past the row's high end is the dropped leading mask of the alignment order
    assert_eq!(run(b"ABCDE", 47, -3, -1, &[frame(-1, 40, 50)]), "ABCDE");
    // a gap column in the row keeps its '-'
    assert_eq!(run(b"AB-DE", 30, 3, 1, &[frame(1, 33, 41)]), "AX-XX");
}
