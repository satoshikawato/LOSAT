use super::super::pattern::AMINOACID_TO_NCBISTDAA;
use super::*;
use crate::config::compat::ScoringMatrix;
use crate::core::blast_util::BLASTAA_SIZE;

fn encode(text: &str) -> Vec<u8> {
    text.bytes()
        .map(|byte| AMINOACID_TO_NCBISTDAA[byte as usize])
        .collect()
}

fn rows_of(query: &[u8]) -> QueryPssm {
    let matrix = crate::utils::matrix::protein_score_rows(ScoringMatrix::Blosum62);
    QueryPssm {
        rows: query
            .iter()
            .map(|&q| {
                let mut row = [0i32; BLASTAA_SIZE];
                for (s, cell) in row.iter_mut().enumerate() {
                    *cell = matrix.score(q, s as u8);
                }
                row
            })
            .collect(),
    }
}

fn costs(x_drop: i32) -> PhiGapCosts {
    PhiGapCosts {
        gap_open: 11,
        gap_extend: 1,
        x_drop,
    }
}

fn consumed(ops: &[GapEditOp]) -> (u32, u32) {
    ops.iter().fold((0, 0), |(q, s), op| match *op {
        GapEditOp::Sub(n) => (q + n, s + n),
        GapEditOp::Del(n) => (q, s + n),
        GapEditOp::Ins(n) => (q + n, s),
    })
}

#[test]
fn self_hit_prelim_and_traceback() {
    let query = encode("MKAGKSTLLW");
    let rows = rows_of(&query);
    let pattern = PhiPattern::compile(b"G-K-S").unwrap();
    let mut scratch = PssmGapScratch::new();
    let prelim = prelim_alignment(&rows, &query, 10, 3, 3, 3, 3, costs(38), &mut scratch);
    // MKA = 5 + 5 + 4 to the left, TLLW = 5 + 4 + 4 + 11 to the right; NCBI's start is one past
    // the extension's.
    assert_eq!(prelim.score, 38);
    assert_eq!((prelim.query_start, prelim.subject_start), (1, 1));
    assert_eq!((prelim.query_stop, prelim.subject_stop), (9, 9));
    let full = traceback(
        &rows,
        &query,
        &query,
        10,
        3,
        3,
        3,
        3,
        &pattern,
        costs(64),
        &mut scratch,
    )
    .unwrap();
    assert_eq!(full.score, 38);
    assert_eq!((full.query_start, full.query_stop), (0, 10));
    assert_eq!(full.edit_script, vec![GapEditOp::Sub(10)]);
}

#[test]
fn prelim_skips_the_left_at_offset_zero() {
    let query = encode("GKSTLLW");
    let rows = rows_of(&query);
    let mut scratch = PssmGapScratch::new();
    let prelim = prelim_alignment(&rows, &query, 7, 0, 0, 3, 3, costs(38), &mut scratch);
    assert_eq!(prelim.score, 24);
    assert_eq!((prelim.query_start, prelim.subject_start), (0, 0));
}

#[test]
fn unequal_wildcard_runs_are_aligned_in_their_band() {
    let pattern = PhiPattern::compile(b"G-x(1,3)-K").unwrap();
    let query = encode("AGLKA");
    let subject = encode("AGLLLKA");
    let rows = rows_of(&query);
    let ops = align_patterns(&rows, &query, &subject, 1, 1, 3, 5, &pattern, 11, 1).unwrap();
    assert_eq!(consumed(&ops), (3, 5));
    // G, then the run of one L against three (one substitution, two deletions; a substitution
    // at either end of the run merges into G's or K's), then K.
    let subs: u32 = ops
        .iter()
        .map(|op| match *op {
            GapEditOp::Sub(n) => n,
            _ => 0,
        })
        .sum();
    assert_eq!(subs, 3);
    assert!(matches!(ops.first(), Some(GapEditOp::Sub(_))));
    assert!(matches!(ops.last(), Some(GapEditOp::Sub(_))));
    // Equal runs are substitutions.
    let ops = align_patterns(&rows, &query, &query, 1, 1, 3, 3, &pattern, 11, 1).unwrap();
    assert_eq!(ops, vec![GapEditOp::Sub(3)]);
}

#[test]
fn banded_runs_at_the_edges() {
    let query = encode("LLLLLL");
    let subject = encode("LLLLLLLLL");
    let rows = rows_of(&query);
    let segments = Segments {
        rows: &rows,
        subject: &subject,
        query_base: 0,
        subject_base: 0,
        gap_open: 11,
        gap_extend: 1,
    };
    let mut ops = Vec::new();
    banded_align(&segments, 0, 3, &mut ops);
    assert_eq!(ops, vec![GapEditOp::Del(3)]);
    let mut ops = Vec::new();
    banded_align(&segments, 2, 0, &mut ops);
    assert_eq!(ops, vec![GapEditOp::Ins(2)]);
    for (start1, start2) in [(2, 5), (5, 2), (4, 9), (6, 1)] {
        let mut ops = Vec::new();
        banded_align(&segments, start1, start2, &mut ops);
        assert_eq!(
            consumed(&ops),
            (start1 as u32, start2 as u32),
            "{start1} {start2}"
        );
    }
}

#[test]
fn very_long_pattern_region() {
    let pattern = PhiPattern::compile(b"[AG]-x(4)-G-K-[ST]-x(0,25)-L").unwrap();
    let query = encode("ALLLLGKSL");
    let subject = encode("ALLLLGKSWWL");
    let rows = rows_of(&query);
    let ops = align_patterns(&rows, &query, &subject, 0, 0, 9, 11, &pattern, 11, 1).unwrap();
    assert_eq!(consumed(&ops), (9, 11));
    assert_eq!(ops.first(), Some(&GapEditOp::Sub(8)));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:497-656
// A long subject's merged HSP aligns the pattern region from one occurrence's start with another
// occurrence's length; where that passes the end of the query or of the subject, NCBI reads past
// the sequence, and LOSAT refuses the search instead.
#[test]
fn a_pattern_region_past_the_end_is_refused() {
    let pattern = PhiPattern::compile(b"G-x(1,3)-K").unwrap();
    let query = encode("AGLKA");
    let subject = encode("AGLLLKA");
    let rows = rows_of(&query);
    for (q_start, s_start, q_len, s_len) in [(1, 1, 5, 5), (1, 3, 3, 5)] {
        let refused = align_patterns(
            &rows, &query, &subject, q_start, s_start, q_len, s_len, &pattern, 11, 1,
        );
        assert!(
            matches!(&refused, Err(PhiError::Refusal(message))
                if message.contains("a pattern region past the end of the query or subject")),
            "{q_start} {s_start} {q_len} {s_len}"
        );
    }
    // The regions that end at the last letters are aligned.
    assert!(align_patterns(&rows, &query, &subject, 1, 1, 4, 6, &pattern, 11, 1).is_ok());
}
