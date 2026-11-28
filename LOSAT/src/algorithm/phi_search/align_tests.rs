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

/// The scores of the gap residue (the zero bytes after the query's letters) against each letter.
fn gap_row() -> [i32; BLASTAA_SIZE] {
    let matrix = crate::utils::matrix::protein_score_rows(ScoringMatrix::Blosum62);
    std::array::from_fn(|s| matrix.score(crate::utils::matrix::ncbistdaa::GAP, s as u8))
}

/// `align_patterns` of the query region `(q_start, q_len)` against the subject region
/// `(s_start, s_len)` with BLOSUM62's rows of `letters` and gap costs 11/1.
fn align(
    letters: &[u8],
    query: &[u8],
    subject: &[u8],
    (q_start, q_len): (usize, usize),
    (s_start, s_len): (usize, usize),
    pattern: &PhiPattern,
) -> Result<Vec<GapEditOp>, PhiError> {
    let rows = rows_of(letters);
    align_patterns(
        &rows,
        &gap_row(),
        query,
        subject,
        q_start,
        s_start,
        q_len,
        s_len,
        pattern,
        11,
        1,
    )
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
        &gap_row(),
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
    let ops = align(&query, &query, &subject, (1, 3), (1, 5), &pattern).unwrap();
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
    let ops = align(&query, &query, &query, (1, 3), (1, 3), &pattern).unwrap();
    assert_eq!(ops, vec![GapEditOp::Sub(3)]);
}

#[test]
fn banded_runs_at_the_edges() {
    let query = encode("LLLLLL");
    let subject = encode("LLLLLLLLL");
    let rows = rows_of(&query);
    let gaps = gap_row();
    let segments = Segments {
        rows: &rows,
        gap_row: &gaps,
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
    let ops = align(&query, &query, &subject, (0, 9), (0, 11), &pattern).unwrap();
    assert_eq!(consumed(&ops), (9, 11));
    assert_eq!(ops.first(), Some(&GapEditOp::Sub(8)));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:497-656
// A long subject's merged HSP aligns the pattern region from one occurrence's start with another
// occurrence's length. Where that passes the end of the query, NCBI reads the two zero bytes after
// its letters (the sentinel and the byte the buffer's allocation adds), and after the subject's
// letters its sentinel; past those it reads outside its buffer (or a database's next sequence),
// and LOSAT refuses the search.
#[test]
fn a_pattern_region_reads_the_zero_bytes_after_the_letters_and_no_further() {
    let pattern = PhiPattern::compile(b"G-x(1,3)-K").unwrap();
    let letters = encode("AGLK");
    let mut query = letters.clone();
    query.extend([0; QUERY_ZEROS_PAST_END]);
    let mut subject = encode("AGLLLK");
    subject.push(0);
    for (q_len, s_len) in [(5, 6), (4, 5), (3, 6)] {
        let aligned = align(&letters, &query, &subject, (1, q_len), (1, s_len), &pattern);
        assert!(aligned.is_ok(), "{q_len} {s_len}: {aligned:?}");
    }
    for (q_len, s_len) in [(6, 5), (4, 7)] {
        let refused = align(&letters, &query, &subject, (1, q_len), (1, s_len), &pattern);
        assert!(
            matches!(&refused, Err(PhiError::Refusal(message))
                if message.contains("a pattern region past the bytes after the query or subject")),
            "{q_len} {s_len}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1579
// The banded aligner scores a query position past the letters with the gap residue's row
// (`gap_row`: NCBI's matrix row of `-`, BLAST_SCORE_MIN against every letter).
#[test]
fn the_banded_aligner_scores_past_the_letters_with_the_gap_row() {
    let rows = rows_of(&encode("AG"));
    let gaps = gap_row();
    assert!(gaps
        .iter()
        .all(|&score| score == crate::utils::matrix::BLAST_SCORE_MIN));
    let subject = encode("LLLK");
    let segments = Segments {
        rows: &rows,
        gap_row: &gaps,
        subject: &subject,
        query_base: 2,
        subject_base: 0,
        gap_open: 11,
        gap_extend: 1,
    };
    assert_eq!(segments.score(1, 1), crate::utils::matrix::BLAST_SCORE_MIN);
    let mut ops = Vec::new();
    banded_align(&segments, 1, 3, &mut ops);
    assert_eq!(consumed(&ops), (1, 3));
}
