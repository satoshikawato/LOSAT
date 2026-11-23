//! Unit tests of the query-anchored section of a PSSM query of the molecule `other` (line PS, PS15b-2
//! fix 1; `write_unreadable_section`): NCBI BLAST+ 2.17.0 merges the HSPs with the query's rows at the
//! width of a nucleotide sequence, and its run stops with the merge's exception where two blocks of
//! one HSP meet on a query row (`X_MERGE_START_ITERATOR`) or the rows order shared segments in two
//! ways (`BUILD_INFINITE_LOOP`). The rows of the real cases are NCBI's own `-outfmt "6 sseqid qstart
//! sstart btop"` rows of psiblast rounds with the molecule `other` PSSMs of the review kit
//! (`/home/kawato/losat-baselines/psiblast-20261025/PS/tmp/rev15b2/inputs/`: `o/tk_other_stdaa.asn`
//! round 1, NCBI's "Start iterator already exists for seq2"; `ps/ps6c_other_iupacaa.asn`
//! `-num_iterations 3 -inclusion_ethresh 10`, rounds 2 (answered) and 3, NCBI's "Infinite loop
//! detected"), with their rows' order and subjects in order of first appearance.

use super::tests::push_btop;
use super::*;

/// The HSPs of `rows` (subject, qstart, sstart, BTOP) of a query whose Bioseq is not a protein.
fn unreadable(rows: &[(u32, u32, u32, &str)]) -> Hsps {
    let mut hsps = Hsps::default();
    for &(subject, qstart, sstart, btop) in rows {
        push_btop(&mut hsps, subject, qstart, sstart, btop);
    }
    hsps.set_query_not_protein();
    hsps
}

/// The section of `hsps`: its bytes, or the refusal's text.
fn section(hsps: &Hsps) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    match write_unreadable_section(&mut out, hsps, "PSIBLAST") {
        Ok(()) => Ok(out),
        Err(error) => {
            assert!(
                error
                    .get_ref()
                    .is_some_and(|inner| inner.is::<UnreadableMergeRefusal>()),
                "{error}"
            );
            Err(error.to_string())
        }
    }
}

#[test]
fn widths_of_a_query_that_is_not_a_protein() {
    let mut protein = Hsps::default();
    assert_eq!(protein.mix_widths(), (1, 1));
    protein.set_query_not_protein();
    assert_eq!(protein.mix_widths(), (3, 1));
    let mut translated = Hsps::of_translated_subjects();
    translated.set_query_not_protein();
    assert_eq!(translated.mix_widths(), (1, 3));
}

#[test]
fn no_hsps_write_nothing() {
    assert_eq!(section(&unreadable(&[])), Ok(Vec::new()));
}

#[test]
fn blocks_on_different_rows_are_merged_and_the_section_is_a_line_feed() {
    // Starts 0, 10 and 11: rows 0, 1 and 2.
    assert_eq!(
        section(&unreadable(&[(0, 1, 1, "10A-10")])),
        Ok(b"\n".to_vec())
    );
    assert_eq!(section(&unreadable(&[(0, 1, 1, "60")])), Ok(b"\n".to_vec()));
}

#[test]
fn two_blocks_of_one_hsp_meeting_on_a_row_are_refused() {
    // Starts 0 and 12 (row 0), spans [0, 30) and [12, 42): the subject gets two starts in one
    // segment.
    let refused = section(&unreadable(&[(0, 1, 1, "10A-A-10")])).unwrap_err();
    assert!(refused.starts_with("-outfmt 1 to 4 with a PSSM query of the molecule other is not supported by LOSAT's PSIBLAST"));
    assert!(!refused.contains("2.17.0+") && !refused.contains("BLAST error:"));
    // The same HSP, its blocks on rows 0 and 1 (start 10): merged.
    assert_eq!(
        section(&unreadable(&[(0, 1, 1, "10A-10")])),
        Ok(b"\n".to_vec())
    );
}

#[test]
fn shared_segments_in_two_orders_are_refused() {
    // HSP 1: rows 0 then 1 (starts 0 and 100); HSP 2: rows 1 then 0 (starts 103 and 114). They
    // share a segment on each row, in the opposite orders of their subjects.
    let query_only = "A-".repeat(50);
    let first = format!("50{query_only}20");
    let hsps = unreadable(&[(0, 1, 1, &first), (1, 104, 1, "10A-10")]);
    assert!(section(&hsps).is_err());
    // Either HSP alone: merged.
    assert_eq!(
        section(&unreadable(&[(0, 1, 1, &first)])),
        Ok(b"\n".to_vec())
    );
    assert_eq!(
        section(&unreadable(&[(1, 104, 1, "10A-10")])),
        Ok(b"\n".to_vec())
    );
}

#[test]
fn the_hsp_that_ncbi_stops_at_in_round_1_of_tk_other_stdaa_is_refused() {
    // BDT63523.1: 18 aligned columns from query start 200, 5 subject letters, 6 aligned columns
    // from query start 218 (row 2 both).
    let hsps = unreadable(&[(0, 201, 484, "1GF1LVVN2EPGD2RQ2KN1TGQE-Y-E-S-M-T1KH1LF2")]);
    assert!(section(&hsps).is_err());
    // A protein query's merge of the same HSP returns.
    let mut protein = Hsps::default();
    push_btop(
        &mut protein,
        0,
        201,
        484,
        "1GF1LVVN2EPGD2RQ2KN1TGQE-Y-E-S-M-T1KH1LF2",
    );
    assert!(merge_groups(&protein).is_ok());
}

/// Round 2 of `ps6c_other_iupacaa.asn -num_iterations 3 -inclusion_ethresh 10` (NCBI answers it).
const IUPAC_ROUND_2: &[(u32, u32, u32, &str)] = &[
    (
        0,
        1,
        4,
        "1MNAI1RKVIEH2TI8TSYTLIKINK3QKEKNIGKGMNQ1HR1LV1VI2SVLKEDTN4GETLGTTEIT3AS1EQVT",
    ),
    (1, 18, 145, "SN1YT1KS2YRQR1ERNLGYGNNF4"),
    (2, 19, 786, "TSYS2NQIL1QKQY1NAGEGN3"),
    (3, 4, 2806, "1RLVLETLF2GS1MR3KR"),
    (
        4,
        22,
        503,
        "1NR2QEQE2GEGE1-V-M-K-KKDHNCLLK1VI1HYSALIEFTLRNYDGT1GSTLGKTKILVFTQ1AQGK1VT",
    ),
    (5, 6, 611, "VIESLC1TI1PLMYFQ1GAKMSD4NK1YRQTQG1"),
];

/// Round 3 of the same run (NCBI stops with "Infinite loop detected").
const IUPAC_ROUND_3: &[(u32, u32, u32, &str)] = &[
    (
        0,
        1,
        4,
        "1MNAI1RKVIEH2TI8TSYTLIKINK3QKEKNIGKGMNQ1HR1LV1VI2SVLKEDTN4GETLGTTEIT3AS1EQVT",
    ),
    (
        1,
        22,
        503,
        "1NR2QEQE2GEGE1-V-M-K-KKDHNCLLK1VI1HYSALIEFTLRNYDGT1GSTLGKTKILVFTQ1AQGK1VT",
    ),
    (
        2,
        23,
        510,
        "NRIV1QEQE3GE1-V-M-K-KKDHNCLLK1VL1HYSALIEFTLRNYDGT1GSTLGKTKILVFTQ1AQGK1VT",
    ),
    (
        3,
        23,
        503,
        "NRIV1QEQE3GE1-V-M-K-KKDHNCLLK1VL1HYSALIEFTLRNYDGT1GSTLGKTKILVFTQ1AQGK1VT",
    ),
    (4, 6, 611, "VIESLC1TI1PLMYFQ1GAKMSD4NK1YRQTQG1"),
    (5, 18, 145, "SN1YT1KS2YRQR1ERNLGYGNNF4"),
    (6, 18, 785, "SETSYS2NQIL1QKQY1NAGEGN3"),
    (7, 4, 2806, "1RLVLETLF2GS1MR3KR"),
    (
        8,
        6,
        278,
        "VIEKLC2GVPLMK2GKKSSNTYYTLI-HKS2YNQN2NLGIGNNQKGHM1LEFNVGKNHV1LTEDTNRIYIGT1GSTIGC1IV",
    ),
    (
        9,
        23,
        591,
        "NR1YRQLQFEY1GQGENE1HN1LGFN1KPHYSVLIEDTDRPYKGSCAGATFGKTD1VLTI1AK",
    ),
];

#[test]
fn the_rounds_of_ps6c_other_iupacaa_are_ncbis() {
    assert_eq!(section(&unreadable(IUPAC_ROUND_2)), Ok(b"\n".to_vec()));
    assert!(section(&unreadable(IUPAC_ROUND_3)).is_err());
}
