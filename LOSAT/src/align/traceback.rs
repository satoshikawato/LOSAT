use super::result::{AlignmentResult, EditOp};

/// Direction for traceback in DP matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TracebackDir {
    /// Diagonal (match/mismatch)
    Diag,
    /// Up (gap in subject / insertion in query)
    Up,
    /// Left (gap in query / deletion from query)
    Left,
    /// Stop (end of alignment)
    Stop,
}

/// Traceback matrix for storing alignment path
pub struct TracebackMatrix {
    data: Vec<TracebackDir>,
    rows: usize,
    cols: usize,
}

impl TracebackMatrix {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![TracebackDir::Stop; rows * cols],
            rows,
            cols,
        }
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> TracebackDir {
        self.data[row * self.cols + col]
    }

    #[inline]
    pub fn set(&mut self, row: usize, col: usize, dir: TracebackDir) {
        self.data[row * self.cols + col] = dir;
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }
}

/// Perform traceback from a given position to reconstruct the alignment
///
/// Returns the edit script and the starting positions (q_start, s_start)
pub fn traceback(
    matrix: &TracebackMatrix,
    end_row: usize,
    end_col: usize,
    query: &[u8],
    subject: &[u8],
) -> (Vec<EditOp>, usize, usize) {
    let mut edit_script = Vec::new();
    let mut row = end_row;
    let mut col = end_col;

    loop {
        let dir = matrix.get(row, col);
        match dir {
            TracebackDir::Diag => {
                if row == 0 || col == 0 {
                    break;
                }
                // Check if match or mismatch
                let q_idx = row - 1;
                let s_idx = col - 1;
                if q_idx < query.len() && s_idx < subject.len() && query[q_idx] == subject[s_idx] {
                    edit_script.push(EditOp::Match);
                } else {
                    edit_script.push(EditOp::Mismatch);
                }
                row -= 1;
                col -= 1;
            }
            TracebackDir::Up => {
                if row == 0 {
                    break;
                }
                edit_script.push(EditOp::Ins);
                row -= 1;
            }
            TracebackDir::Left => {
                if col == 0 {
                    break;
                }
                edit_script.push(EditOp::Del);
                col -= 1;
            }
            TracebackDir::Stop => {
                break;
            }
        }
    }

    // Reverse to get forward order
    edit_script.reverse();

    (edit_script, row, col)
}

/// Compute alignment statistics from two sequences and their alignment
pub fn compute_alignment_stats(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
) -> (usize, usize, usize, usize) {
    // Simple ungapped alignment statistics
    let q_len = q_end.saturating_sub(q_start);
    let s_len = s_end.saturating_sub(s_start);
    let alignment_len = q_len.max(s_len);

    let mut matches = 0;
    let mut mismatches = 0;

    let min_len = q_len.min(s_len);
    for i in 0..min_len {
        let q_idx = q_start + i;
        let s_idx = s_start + i;
        if q_idx < query.len() && s_idx < subject.len() {
            if query[q_idx] == subject[s_idx] {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }

    // Remaining positions are gaps
    let gap_opens = if q_len != s_len { 1 } else { 0 };

    (matches, mismatches, gap_opens, alignment_len)
}

/// Compute statistics from an edit script
pub fn stats_from_edit_script(edit_script: &[EditOp]) -> (usize, usize, usize, usize) {
    let mut matches = 0;
    let mut mismatches = 0;
    let mut gap_opens = 0;
    let mut prev_op: Option<EditOp> = None;

    for &op in edit_script {
        match op {
            EditOp::Match => matches += 1,
            EditOp::Mismatch => mismatches += 1,
            EditOp::Ins => {
                if prev_op != Some(EditOp::Ins) {
                    gap_opens += 1;
                }
            }
            EditOp::Del => {
                if prev_op != Some(EditOp::Del) {
                    gap_opens += 1;
                }
            }
        }
        prev_op = Some(op);
    }

    (matches, mismatches, gap_opens, edit_script.len())
}

/// Generate edit script from two aligned sequences (simple ungapped case)
pub fn generate_ungapped_edit_script(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    length: usize,
) -> Vec<EditOp> {
    let mut script = Vec::with_capacity(length);

    for i in 0..length {
        let q_idx = q_start + i;
        let s_idx = s_start + i;

        if q_idx < query.len() && s_idx < subject.len() {
            if query[q_idx] == subject[s_idx] {
                script.push(EditOp::Match);
            } else {
                script.push(EditOp::Mismatch);
            }
        }
    }

    script
}

/// Create an AlignmentResult from traceback
pub fn alignment_from_traceback(
    matrix: &TracebackMatrix,
    end_row: usize,
    end_col: usize,
    query: &[u8],
    subject: &[u8],
    score: i32,
) -> AlignmentResult {
    let (edit_script, start_row, start_col) = traceback(matrix, end_row, end_col, query, subject);

    AlignmentResult::with_traceback(
        start_row + 1, // Convert to 1-based
        end_row,       // Already 1-based (row index = position)
        start_col + 1, // Convert to 1-based
        end_col,       // Already 1-based
        score,
        edit_script,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_traceback_matrix() {
        let mut matrix = TracebackMatrix::new(5, 5);
        matrix.set(2, 3, TracebackDir::Diag);
        assert_eq!(matrix.get(2, 3), TracebackDir::Diag);
        assert_eq!(matrix.get(0, 0), TracebackDir::Stop);
    }

    #[test]
    fn test_stats_from_edit_script() {
        let script = vec![
            EditOp::Match,
            EditOp::Match,
            EditOp::Mismatch,
            EditOp::Ins,
            EditOp::Ins,
            EditOp::Match,
        ];

        let (matches, mismatches, gap_opens, len) = stats_from_edit_script(&script);
        assert_eq!(matches, 3);
        assert_eq!(mismatches, 1);
        assert_eq!(gap_opens, 1);
        assert_eq!(len, 6);
    }

    #[test]
    fn test_ungapped_edit_script() {
        let query = b"ACGTACGT";
        let subject = b"ACGTTCGT";

        let script = generate_ungapped_edit_script(query, subject, 0, 0, 8);

        assert_eq!(script.len(), 8);
        assert_eq!(script[0], EditOp::Match); // A=A
        assert_eq!(script[4], EditOp::Mismatch); // A!=T
    }
}
