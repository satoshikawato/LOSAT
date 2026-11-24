//! Unit tests of the query-anchored section of a PSSM query of the molecule `other` (line PS,
//! PS15b-2 fix 1; `write_unreadable_section`): NCBI BLAST+ 2.17.0 merges the HSPs with the query's
//! rows at the width of a nucleotide sequence, and its run stops where the merge fails: with an
//! exception where two blocks of one subject row meet on a query row (`X_MERGE_START_ITERATOR`),
//! where the rows order shared segments in two ways (`BUILD_INFINITE_LOOP`), where a later pass
//! finds no row of its HSP (`X_MERGE_PROPER_ROW`) or where `Build` reads past the end of a row's
//! starts (`MixError::NcbiReadsPastEnd`), with SIGSEGV where a query letter against a gap of the
//! subject opens a row of its own (`MixError::NcbiNullSecondRow`), and with a read of freed memory
//! where `FillUnalignedRegions` replaces the last reference to a listed segment
//! (`MixError::NcbiFreesSegment`); a length below zero (`MixError::NcbiWrapsLength`) is refused. Each expectation is NCBI's own run: the rows are NCBI's `-outfmt
//! "6 sseqid qstart sstart btop"` rows (from the PSSM's molecule `aa` twin, whose search is the
//! same), in their order, subjects numbered in order of first appearance, of psiblast with the
//! review kit's molecule `other` PSSMs
//! (`/home/kawato/losat-baselines/psiblast-20261025/PS/tmp/rev15b2/inputs/`:
//! `o/tk_other_stdaa.asn`, `ps/ps6c_other_iupacaa.asn`) and the subjects of
//! `/home/kawato/losat-baselines/psiblast-20261025/PS/tmp/fix15b2_t0/` (the PSSM's query with
//! letters deleted or inserted, and pieces of it), whose `-outfmt 1` NCBI answered (exit status 0),
//! ended with the merge's exception (255) or with SIGSEGV (139).

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

/// The section of `hsps`: its bytes, or how NCBI's merge fails (the refusal's), before which the
/// section writes nothing.
fn section(hsps: &Hsps) -> Result<Vec<u8>, MergeFailure> {
    let mut out = Vec::new();
    match write_unreadable_section(&mut out, hsps, "PSIBLAST") {
        Ok(()) => Ok(out),
        Err(error) => {
            let Some(refusal) = error
                .get_ref()
                .and_then(|inner| inner.downcast_ref::<UnreadableMergeRefusal>())
            else {
                panic!("{error}");
            };
            assert!(out.is_empty());
            Err(refusal.failure)
        }
    }
}

/// The section of `rows` (`unreadable`, `section`).
fn rows_section(rows: &[(u32, u32, u32, &str)]) -> Result<Vec<u8>, MergeFailure> {
    section(&unreadable(rows))
}

/// The section of a merge that returns: a line feed.
fn line_feed() -> Result<Vec<u8>, MergeFailure> {
    Ok(b"\n".to_vec())
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
    assert_eq!(rows_section(&[]), Ok(Vec::new()));
}

#[test]
fn blocks_on_different_rows_are_merged_and_the_section_is_a_line_feed() {
    assert_eq!(rows_section(&[(0, 1, 1, "60")]), line_feed());
    // Starts 0 and 10: rows 0 and 1; the subject's letter alone goes to the subject's pass.
    assert_eq!(rows_section(&[(0, 1, 1, "10-A10")]), line_feed());
    // NCBI: `self`, `ins100`, `ins101` (a letter inserted in the subject).
    assert_eq!(rows_section(&[(0, 1, 1, "398")]), line_feed());
    assert_eq!(rows_section(&[(0, 1, 1, "100-W298")]), line_feed());
    assert_eq!(rows_section(&[(0, 1, 1, "101-W297")]), line_feed());
}

#[test]
fn two_blocks_of_one_hsp_meeting_on_a_row_are_ncbis_exception() {
    // Starts 0 and 12 (row 0), spans [0, 30) and [12, 42): the subject gets two starts in one
    // segment.
    assert_eq!(
        rows_section(&[(0, 1, 1, "10A-A-10")]),
        Err(MergeFailure::Exception)
    );
    // NCBI: `del101` (starts 0 and 102), `del2gap` (0 and 102), `del100_del200` (0 and 201, after
    // a query letter alone on row 1, which the matches of two rows come before).
    for btop in ["101T-296", "100F-T-296", "100F-99K-197"] {
        assert_eq!(
            rows_section(&[(0, 1, 1, btop)]),
            Err(MergeFailure::Exception),
            "{btop}"
        );
    }
}

#[test]
fn a_query_letter_against_a_gap_on_a_row_of_its_own_is_ncbis_crash() {
    // Starts 0 and 11 (rows 0 and 2) and the query's letter alone at 10 (row 1).
    assert_eq!(
        rows_section(&[(0, 1, 1, "10A-10")]),
        Err(MergeFailure::Crash)
    );
    // NCBI: `del100`, `del103` (rows 0, 1 alone, 2), `two_self` and `two_ins` (another subject's
    // blocks on rows 0 and 2).
    assert_eq!(
        rows_section(&[(0, 1, 1, "100F-297")]),
        Err(MergeFailure::Crash)
    );
    assert_eq!(
        rows_section(&[(0, 1, 1, "103L-294")]),
        Err(MergeFailure::Crash)
    );
    assert_eq!(
        rows_section(&[(0, 1, 1, "398"), (1, 1, 1, "100F-297")]),
        Err(MergeFailure::Crash)
    );
    assert_eq!(
        rows_section(&[(0, 1, 1, "101-W297"), (1, 1, 1, "100F-297")]),
        Err(MergeFailure::Crash)
    );
}

#[test]
fn a_query_letter_against_a_gap_on_an_open_row_is_merged() {
    // NCBI: `del102` (the letter alone at 102 in row 0, inside the first block's span), `del3gap`
    // (at 100 in row 1, the row of the block after it).
    assert_eq!(rows_section(&[(0, 1, 1, "102D-295")]), line_feed());
    assert_eq!(rows_section(&[(0, 1, 1, "100F-T-D-295")]), line_feed());
    // NCBI: `two_cover` and `two_frag2`: another subject's block opened row 1.
    assert_eq!(
        rows_section(&[
            (0, 1, 1, "100F-297"),
            (1, 2, 1, "299"),
            (1, 317, 95, "3LQPFQFPT4LV")
        ]),
        line_feed()
    );
    assert_eq!(
        rows_section(&[
            (0, 1, 1, "100F-297"),
            (1, 3, 1, "298"),
            (1, 317, 94, "3LQPFQFPT4LV")
        ]),
        line_feed()
    );
}

/// The HSPs of subject `s1` of `t2a.faa` (its query's first 121 letters with a letter inserted
/// after 61).
const T2A_S1: &[(u32, u32, u32, &str)] = &[
    (0, 1, 1, "61-W60"),
    (0, 195, 108, "IV1KRKLKAND3LRVI1LI"),
    (0, 317, 97, "3LQPFQFPT4LV"),
    (0, 323, 88, "1DELV1LDLVMI1LI1"),
];

#[test]
fn shared_segments_in_two_orders_are_ncbis_exception() {
    // HSP 1: row 0 then row 1 (starts 0 and 10); HSP 2: row 1 then row 0 (starts 13 and 18). They
    // share a segment on each row, in the opposite orders of their subjects.
    assert_eq!(
        rows_section(&[(0, 1, 1, "10-A10"), (1, 14, 1, "5-A10")]),
        Err(MergeFailure::Exception)
    );
    assert_eq!(rows_section(&[(1, 14, 1, "5-A10")]), line_feed());
    // NCBI: `t2a -evalue 1e-10` (the same with blocks of 61, 60, 32 and 60 letters; Infinite loop
    // detected), and `t2none -evalue 1e-10`, whose second HSP's row 0 block starts after the first
    // HSP's row 0 span (no shared segment there).
    assert_eq!(
        rows_section(&[(0, 1, 1, "61-W60"), (1, 101, 1, "32-W60")]),
        Err(MergeFailure::Exception)
    );
    assert_eq!(
        rows_section(&[(0, 1, 1, "61-W60"), (1, 182, 1, "32-W60")]),
        line_feed()
    );
}

#[test]
fn a_later_pass_without_the_row_of_its_hsp_is_ncbis_exception() {
    // NCBI: `t2a_s1` (Proper row not found for the match): the first match of row 2 (start 316, the
    // third HSP's) gives the subject's first row to that HSP, so the letter alone of the first HSP
    // finds no row of its own; `t2none` (the same, with another subject).
    assert_eq!(rows_section(T2A_S1), Err(MergeFailure::Exception));
    let mut t2none = T2A_S1.to_vec();
    t2none.push((1, 182, 1, "32-W60"));
    assert_eq!(rows_section(&t2none), Err(MergeFailure::Exception));
    // NCBI: `t2a_s2` and `t2a` (Start iterator already exists for seq2): the first HSP's row taken
    // by a later HSP of the subject, whose blocks then meet its own.
    let t2a_s2 = [
        (0, 101, 1, "32-W60"),
        (0, 195, 7, "IV1KRKLKAND3LRVI1LIEAGAGL1R-C-1KT1TD1"),
    ];
    assert_eq!(rows_section(&t2a_s2), Err(MergeFailure::Exception));
    let mut t2a = T2A_S1.to_vec();
    t2a.extend(t2a_s2.iter().map(|&(_, q, s, btop)| (1, q, s, btop)));
    assert_eq!(rows_section(&t2a), Err(MergeFailure::Exception));
}

#[test]
fn the_hsp_that_ncbi_stops_at_in_round_1_of_tk_other_stdaa_is_refused() {
    // BDT63523.1: 18 aligned columns from query start 200, 5 subject letters, 6 aligned columns
    // from query start 218 (row 2 both).
    let hsps = unreadable(&[(0, 201, 484, "1GF1LVVN2EPGD2RQ2KN1TGQE-Y-E-S-M-T1KH1LF2")]);
    assert_eq!(section(&hsps), Err(MergeFailure::Exception));
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
    assert_eq!(rows_section(IUPAC_ROUND_2), line_feed());
    assert_eq!(rows_section(IUPAC_ROUND_3), Err(MergeFailure::Exception));
}

/// Round 1 of `tk_other_stdaa.asn -subject s.faa -num_iterations 2` of the differential corpus's
/// case z147 (`/home/kawato/losat-baselines/psiblast-20261025/PS/tmp/fix15b2_dz/cases_z/z147`,
/// letters 195-313 of the PSSM's query; row `db.ps15b2.p.othr.past.1`), which NCBI answers.
const Z147_ROUND_1: &[(u32, u32, u32, &str)] = &[
    (0, 195, 1, "119"),
    (0, 11, 6, "MNAK1RKVLEVLVVLTE1PGMDFRAC4YQLA1NLILYLQTQNEK1"),
    (0, 122, 65, "2GDTHSADAQHQLML2PAI-H-KNLRLWPEYVTCNSSK1VKKQ1"),
];

/// Round 2 of the same run: NCBI stops in `CAlnMixSegments::Build` with "CObject is corrupted"
/// (exit status 255, 3 runs of 3). The first match of the frame-1 row (the second HSP's first
/// block) gives the subject's first row to the second HSP, whose last block then lands on that
/// row inside the first HSP's segment.
const Z147_ROUND_2: &[(u32, u32, u32, &str)] = &[
    (0, 195, 1, "119"),
    (0, 107, 1, "VI1RKLKAKDN3RLIV1ILAEAGLG1-R-C1TK1DT1"),
];

#[test]
fn a_row_walked_past_its_end_is_ncbis_exception() {
    assert_eq!(rows_section(Z147_ROUND_1), line_feed());
    let hsps = unreadable(Z147_ROUND_2);
    let merged = merge_groups(&hsps);
    assert!(
        matches!(merged, Err(aln_mix::MixError::NcbiReadsPastEnd)),
        "{:?}",
        merged.err()
    );
    assert_eq!(section(&hsps), Err(MergeFailure::Exception));
}

/// xorshift64*: the cases of `random_sections_end_as_ncbis_merge_can`, the same in every run.
struct Cases(u64);

impl Cases {
    fn below(&mut self, n: u32) -> u32 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        ((x.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 33) % u64::from(n)) as u32
    }

    /// A BTOP of 1 to 4 aligned blocks of 1 to 40 columns, with 1 to 4 letters of the query or of
    /// the subject alone between two blocks.
    fn btop(&mut self) -> String {
        let mut btop = String::new();
        for block in 0..1 + self.below(4) {
            if block > 0 {
                let gap = if self.below(2) == 0 { "A-" } else { "-A" };
                for _ in 0..1 + self.below(4) {
                    btop.push_str(gap);
                }
            }
            btop.push_str(&(1 + self.below(40)).to_string());
        }
        btop
    }
}

/// Case 586 of `random_sections_end_as_ncbis_merge_can`: a fill of the frame-1 row starts on a
/// start of that row, whose segment it frees while `m_Segments` lists it (no NCBI run: a read of
/// freed memory has no defined result).
#[test]
fn a_fill_that_frees_a_listed_segment_is_refused() {
    let hsps = unreadable(&[
        (0, 31, 134, "12"),
        (0, 42, 44, "8-A-A1A-A-A-A-8"),
        (0, 88, 96, "40"),
        (1, 39, 109, "16"),
    ]);
    let merged = merge_groups(&hsps);
    assert!(
        matches!(merged, Err(aln_mix::MixError::NcbiFreesSegment)),
        "{:?}",
        merged.err()
    );
    assert_eq!(section(&hsps), Err(MergeFailure::Freed));
}

/// Case 2890 of `random_sections_end_as_ncbis_merge_can`: a match's walk on its query row finds
/// the next start before its own position, where NCBI's length of the new segment wraps (no NCBI
/// run: NCBI then runs forever or goes on with lengths near 2^32).
#[test]
fn a_length_below_zero_is_refused() {
    let hsps = unreadable(&[
        (0, 115, 77, "36A-3"),
        (0, 120, 69, "10-A-A-A29A-18A-A-A-A-29"),
        (1, 116, 150, "23"),
        (1, 8, 4, "9"),
        (1, 65, 4, "38-A-A32"),
        (2, 46, 1, "4A-A-4A-18-A-A-A-A12"),
        (2, 116, 109, "9-A-A-A-A11"),
        (3, 102, 126, "24-A-A-A-A33A-A-A-13"),
        (3, 4, 67, "39"),
        (3, 116, 89, "15A-35A-5A-A-A-A-39"),
    ]);
    let merged = merge_groups(&hsps);
    assert!(
        matches!(merged, Err(aln_mix::MixError::NcbiWrapsLength)),
        "{:?}",
        merged.err()
    );
    assert_eq!(section(&hsps), Err(MergeFailure::Wraps));
}

#[test]
fn random_sections_end_as_ncbis_merge_can() {
    // Each merge of the HSPs that a search can give such a query returns (a line feed) or ends as
    // NCBI's does: never at a LOSAT internal error or a panic. 1 to 4 subjects of 1 to 3 HSPs
    // each, anywhere on the first 150 letters of the query and the subject.
    let mut cases = Cases(0x5eed_0015_b2f1_0001);
    let mut ends = [0usize; 5];
    let mut unexpected = Vec::new();
    for case in 0..40_000 {
        let mut hsps = Vec::new();
        for subject in 0..1 + cases.below(4) {
            for _ in 0..1 + cases.below(3) {
                let (qstart, sstart) = (1 + cases.below(150), 1 + cases.below(150));
                hsps.push((subject, qstart, sstart, cases.btop()));
            }
        }
        let rows: Vec<(u32, u32, u32, &str)> = hsps
            .iter()
            .map(|(subject, qstart, sstart, btop)| (*subject, *qstart, *sstart, btop.as_str()))
            .collect();
        let mut out = Vec::new();
        let written = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            write_unreadable_section(&mut out, &unreadable(&rows), "PSIBLAST")
        }));
        let end = match written {
            Ok(Ok(())) => {
                assert_eq!(out, b"\n", "case {case}: {rows:?}");
                0
            }
            Ok(Err(error)) => match error
                .get_ref()
                .and_then(|inner| inner.downcast_ref::<UnreadableMergeRefusal>())
                .map(|refusal| refusal.failure)
            {
                Some(MergeFailure::Exception) => 1,
                Some(MergeFailure::Crash) => 2,
                Some(MergeFailure::Freed) => 3,
                Some(MergeFailure::Wraps) => 4,
                None => {
                    unexpected.push(format!("case {case}: {rows:?}: {error}"));
                    continue;
                }
            },
            Err(_) => {
                unexpected.push(format!("case {case}: {rows:?}: a panic"));
                continue;
            }
        };
        ends[end] += 1;
    }
    // Every case ends as one of NCBI's ends, and each end is met (line feed, exception, crash,
    // freed memory, a length below zero).
    assert!(
        unexpected.is_empty(),
        "{} cases: {:?}",
        unexpected.len(),
        &unexpected[..unexpected.len().min(8)]
    );
    assert!(ends.iter().all(|&count| count > 0), "{ends:?}");
}
