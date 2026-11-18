//! Unit tests of `output_seq`'s mask spans under `-html` (line RG, item T4a0; reachable from item
//! T4a): NCBI BLAST+ 2.17.0's lines of the probes of survey T4 §1.5
//! (`/home/kawato/losat-baselines/remaining-gaps-20261112/T/t4/out/<id>.h.out`), and the cases that
//! `x_OutputSeq`'s loop gives where no probe reaches (the range `(0, 0)` that NCBI drops,
//! overlapping masks).

use super::*;

/// One line of `sequence` (columns `start..start + len`) with masks `masks` (columns, frame 0) on a
/// row of frame 0.
fn line(sequence: &[u8], start: usize, len: usize, masks: &[(i64, i64)], html: bool) -> String {
    let masks: Vec<(i64, i64, i8)> = masks.iter().map(|&(from, to)| (from, to, 0)).collect();
    let mut out = Vec::new();
    output_seq(&mut out, sequence, start, len, Some((&masks, 0)), html);
    String::from_utf8(out).unwrap()
}

const OPEN: &str = "<font color=\"#000000\">";
const CLOSE: &str = "</font>";

/// `bn.lc.ll30.f0.h` (`-lcase_masking`, `-line_length 30`): a subject mask of columns 20 to 49
/// spans two lines; on the first it runs to the line's end, on the second it closes at its end.
#[test]
fn a_mask_over_two_lines_opens_on_the_first_and_closes_on_the_second() {
    let row = b"GGTGTACTTTTCTCCCCCTCGTTGTAGCCTAAGGGTTTAACCCGCACGCCGTGGCTTGAG";
    let masks = [(20, 49)];
    assert_eq!(
        line(row, 0, 30, &masks, true),
        format!("GGTGTACTTTTCTCCCCCTC{OPEN}gttgtagcct{CLOSE}")
    );
    assert_eq!(
        line(row, 30, 30, &masks, true),
        format!("{OPEN}aagggtttaacccgcacgcc{CLOSE}GTGGCTTGAG")
    );
    // Without `-html` the same letters are in lower case, without spans.
    assert_eq!(
        line(row, 0, 30, &masks, false),
        "GGTGTACTTTTCTCCCCCTCgttgtagcct"
    );
    assert_eq!(
        line(row, 30, 30, &masks, false),
        "aagggtttaacccgcacgccGTGGCTTGAG"
    );
}

/// `bn.dust.ll7.f0.h` (DUST, `-line_length 7`): `Query  196  AGCAT<font ...>aa</font>  202` and
/// `Query  231  <font ...>aCGTCTG</font>  237`: a mask that ends at a line's first column opens
/// there and runs to the line's end (`to > 0` is false).
#[test]
fn a_mask_that_ends_at_a_lines_first_column_runs_to_its_end() {
    let row = b"AGCATAAAAAAAAACGTCTG";
    let masks = [(5, 13)];
    assert_eq!(
        line(row, 0, 7, &masks, true),
        format!("AGCAT{OPEN}aa{CLOSE}")
    );
    assert_eq!(
        line(row, 7, 7, &masks, true),
        format!("{OPEN}aaaaaaa{CLOSE}")
    );
    let row = b"AAAAAAAACGTCTG";
    assert_eq!(
        line(row, 7, 7, &[(0, 7)], true),
        format!("{OPEN}aCGTCTG{CLOSE}")
    );
}

/// `bn.lc.ll30.f0.h`: `Sbjct  151  -<font ...>cgcg...</font>  179`: a gap before the mask stays
/// outside its span.
#[test]
fn a_gap_before_a_mask_stays_outside_its_span() {
    let row = b"-CGCGTCTGGCGTACGATGACTAACGCTCT";
    assert_eq!(
        line(row, 0, 30, &[(1, 29)], true),
        format!("-{OPEN}cgcgtctggcgtacgatgactaacgctct{CLOSE}")
    );
}

/// A mask from the alignment's first column that does not end on the first line, or covers that
/// column alone, is the range `(0, 0)`, which NCBI drops: its letters are in lower case without a
/// span (survey T4 §1.5; from the source, no probe reaches it).
#[test]
fn a_mask_from_the_first_column_has_no_span_on_the_first_line() {
    let row = b"ACGTCTGACGTCTG";
    assert_eq!(line(row, 0, 7, &[(0, 9)], true), "acgtctg");
    assert_eq!(
        line(row, 7, 7, &[(0, 9)], true),
        format!("{OPEN}acg{CLOSE}TCTG")
    );
    assert_eq!(line(row, 0, 7, &[(0, 0)], true), "aCGTCTG");
    // A mask from column 0 that ends on the first line has its span.
    assert_eq!(
        line(row, 0, 7, &[(0, 2)], true),
        format!("{OPEN}acg{CLOSE}TCTG")
    );
}

/// Two masks on a line each get a span; overlapping masks close at the first end that comes after
/// an opening, and the rest of the later mask is in lower case without a span (NCBI's loop).
#[test]
fn masks_on_one_line_follow_ncbis_loop() {
    let row = b"ABCDEFGHIJ";
    assert_eq!(
        line(row, 0, 10, &[(2, 3), (6, 8)], true),
        format!("AB{OPEN}cd{CLOSE}EF{OPEN}ghi{CLOSE}J")
    );
    assert_eq!(
        line(row, 0, 10, &[(2, 5), (4, 7)], true),
        format!("AB{OPEN}cdef{CLOSE}ghIJ")
    );
    // A mask of another frame is not the row's.
    let masks = [(2, 3, 1)];
    let mut out = Vec::new();
    output_seq(&mut out, row, 0, 10, Some((&masks, 0)), true);
    assert_eq!(out, b"ABCDEFGHIJ");
}
