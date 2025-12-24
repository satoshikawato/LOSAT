/// Result of a sequence alignment with full statistics
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Query start position (1-based)
    pub q_start: usize,
    /// Query end position (1-based, inclusive)
    pub q_end: usize,
    /// Subject start position (1-based)
    pub s_start: usize,
    /// Subject end position (1-based, inclusive)
    pub s_end: usize,
    /// Raw alignment score
    pub score: i32,
    /// Number of identical positions
    pub matches: usize,
    /// Number of mismatched positions
    pub mismatches: usize,
    /// Number of gap openings
    pub gap_opens: usize,
    /// Total alignment length (number of columns including gaps)
    pub alignment_len: usize,
    /// Optional edit script for traceback
    pub edit_script: Option<Vec<EditOp>>,
}

impl AlignmentResult {
    /// Create a new alignment result
    pub fn new(
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        score: i32,
        matches: usize,
        mismatches: usize,
        gap_opens: usize,
        alignment_len: usize,
    ) -> Self {
        Self {
            q_start,
            q_end,
            s_start,
            s_end,
            score,
            matches,
            mismatches,
            gap_opens,
            alignment_len,
            edit_script: None,
        }
    }

    /// Create alignment result with edit script
    pub fn with_traceback(
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        score: i32,
        edit_script: Vec<EditOp>,
    ) -> Self {
        let stats = compute_stats_from_edit_script(&edit_script);
        Self {
            q_start,
            q_end,
            s_start,
            s_end,
            score,
            matches: stats.matches,
            mismatches: stats.mismatches,
            gap_opens: stats.gap_opens,
            alignment_len: stats.alignment_len,
            edit_script: Some(edit_script),
        }
    }

    /// Calculate percent identity
    pub fn identity(&self) -> f64 {
        if self.alignment_len == 0 {
            return 0.0;
        }
        100.0 * (self.matches as f64) / (self.alignment_len as f64)
    }

    /// Calculate query coverage
    pub fn query_coverage(&self) -> f64 {
        let q_len = self.q_end.saturating_sub(self.q_start) + 1;
        if q_len == 0 {
            return 0.0;
        }
        100.0 * (self.alignment_len as f64) / (q_len as f64)
    }

    /// Get the number of gaps (total gap positions, not gap openings)
    pub fn gaps(&self) -> usize {
        self.alignment_len
            .saturating_sub(self.matches + self.mismatches)
    }
}

/// Edit operation for traceback
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditOp {
    /// Match (identical residues)
    Match,
    /// Mismatch (different residues)
    Mismatch,
    /// Insertion in query (gap in subject)
    Ins,
    /// Deletion from query (gap in query)
    Del,
}

/// Statistics computed from edit script
struct EditStats {
    matches: usize,
    mismatches: usize,
    gap_opens: usize,
    alignment_len: usize,
}

/// Compute alignment statistics from edit script
fn compute_stats_from_edit_script(edit_script: &[EditOp]) -> EditStats {
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

    EditStats {
        matches,
        mismatches,
        gap_opens,
        alignment_len: edit_script.len(),
    }
}

/// Builder for creating alignment results incrementally
#[derive(Debug, Default)]
pub struct AlignmentBuilder {
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    score: i32,
    matches: usize,
    mismatches: usize,
    gap_opens: usize,
    alignment_len: usize,
    edit_script: Vec<EditOp>,
}

impl AlignmentBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_positions(
        mut self,
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
    ) -> Self {
        self.q_start = q_start;
        self.q_end = q_end;
        self.s_start = s_start;
        self.s_end = s_end;
        self
    }

    pub fn set_score(mut self, score: i32) -> Self {
        self.score = score;
        self
    }

    pub fn add_match(mut self) -> Self {
        self.matches += 1;
        self.alignment_len += 1;
        self.edit_script.push(EditOp::Match);
        self
    }

    pub fn add_mismatch(mut self) -> Self {
        self.mismatches += 1;
        self.alignment_len += 1;
        self.edit_script.push(EditOp::Mismatch);
        self
    }

    pub fn add_insertion(mut self) -> Self {
        if self.edit_script.last() != Some(&EditOp::Ins) {
            self.gap_opens += 1;
        }
        self.alignment_len += 1;
        self.edit_script.push(EditOp::Ins);
        self
    }

    pub fn add_deletion(mut self) -> Self {
        if self.edit_script.last() != Some(&EditOp::Del) {
            self.gap_opens += 1;
        }
        self.alignment_len += 1;
        self.edit_script.push(EditOp::Del);
        self
    }

    pub fn build(self) -> AlignmentResult {
        AlignmentResult {
            q_start: self.q_start,
            q_end: self.q_end,
            s_start: self.s_start,
            s_end: self.s_end,
            score: self.score,
            matches: self.matches,
            mismatches: self.mismatches,
            gap_opens: self.gap_opens,
            alignment_len: self.alignment_len,
            edit_script: if self.edit_script.is_empty() {
                None
            } else {
                Some(self.edit_script)
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignment_result_identity() {
        let result = AlignmentResult::new(1, 100, 1, 100, 200, 90, 10, 0, 100);
        assert!((result.identity() - 90.0).abs() < 0.001);
    }

    #[test]
    fn test_edit_script_stats() {
        let script = vec![
            EditOp::Match,
            EditOp::Match,
            EditOp::Mismatch,
            EditOp::Ins,
            EditOp::Ins,
            EditOp::Match,
            EditOp::Del,
        ];

        let stats = compute_stats_from_edit_script(&script);
        assert_eq!(stats.matches, 3);
        assert_eq!(stats.mismatches, 1);
        assert_eq!(stats.gap_opens, 2); // One insertion run, one deletion
        assert_eq!(stats.alignment_len, 7);
    }

    #[test]
    fn test_alignment_builder() {
        let result = AlignmentBuilder::new()
            .set_positions(1, 10, 1, 10)
            .set_score(50)
            .add_match()
            .add_match()
            .add_mismatch()
            .add_insertion()
            .add_insertion()
            .add_match()
            .build();

        assert_eq!(result.matches, 3);
        assert_eq!(result.mismatches, 1);
        assert_eq!(result.gap_opens, 1);
        assert_eq!(result.alignment_len, 6);
    }
}
