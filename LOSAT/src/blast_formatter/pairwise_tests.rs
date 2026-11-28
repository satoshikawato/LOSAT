use super::*;

// A stored gap residue is a `-` at a column where the row has letters; a gap column's `-` is
// not one.
#[test]
fn stored_gap_columns_are_letters_shown_as_gaps() {
    let ds = PairDenseg {
        starts: vec![0, 0, -1, 3, 3, 5],
        lens: vec![3, 2, 2],
        strands: [Strand::Unknown, Strand::Unknown],
        first_minus: false,
    };
    // Query row: letters, then a gap of 2, then letters; its `-` at column 1 is a residue.
    assert_eq!(stored_gap_columns(&ds, b"A-C--DE", 0), vec![1]);
    // Subject row: letters everywhere; its `-` at columns 3 and 6 are residues.
    assert_eq!(stored_gap_columns(&ds, b"ACD-EF-", 1), vec![3, 6]);
    assert!(stored_gap_columns(&ds, b"ACDEFGH", 1).is_empty());
}

// The subject letters that the rows read: from the first to the last letter of row 1.
#[test]
fn subject_range_spans_row_one() {
    let ds = PairDenseg {
        starts: vec![10, 20, 13, -1, 15, 23],
        lens: vec![3, 2, 4],
        strands: [Strand::Plus, Strand::Plus],
        first_minus: false,
    };
    assert_eq!(subject_range(&ds), Some((20, 26)));
}

/// A sink whose writes fail from the `fail_at`th on, recording each failure the adapter routes to
/// `alignments_failed` (`flushing`).
struct FailingSink {
    writes: usize,
    fail_at: usize,
    routed: Vec<bool>,
}

impl ReportSink for FailingSink {
    fn write(&mut self, _bytes: &[u8]) -> Result<()> {
        self.writes += 1;
        if self.writes >= self.fail_at {
            return Err(anyhow::anyhow!("write failed"));
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        self.write(b"")
    }

    fn post(&mut self, _lines: &[u8]) {}

    fn badbit(&mut self) -> anyhow::Error {
        anyhow::anyhow!("badbit")
    }

    fn alignments_failed(&mut self, error: anyhow::Error, flushing: bool) -> anyhow::Error {
        self.routed.push(flushing);
        error
    }
}

#[test]
fn a_failure_in_the_alignments_ends_as_there_and_elsewhere_as_the_report() {
    use crate::report::ostream::DisplaySeqalign;
    // Before the alignments (the description table): the report's end, nothing routed.
    let mut sink = FailingSink {
        writes: 0,
        fail_at: 1,
        routed: Vec::new(),
    };
    {
        let cell = SinkCell {
            sink: RefCell::new(&mut sink),
            failed: RefCell::new(None),
            alignments: Cell::new(false),
        };
        assert!(StreamOut(&cell).write_all(b"descriptions\n").is_err());
        assert!(cell.failed.borrow().is_some());
    }
    assert!(sink.routed.is_empty());
    // An insertion of an alignment, then (another run) the flush at the alignments' end.
    for (fail_at, routed) in [(2, vec![false]), (3, vec![true])] {
        let mut sink = FailingSink {
            writes: 0,
            fail_at,
            routed: Vec::new(),
        };
        {
            let cell = SinkCell {
                sink: RefCell::new(&mut sink),
                failed: RefCell::new(None),
                alignments: Cell::new(false),
            };
            let mut out = StreamOut(&cell);
            let mut display = DisplaySeqalign::default();
            out.write_all(b"heading\n").unwrap();
            display.enter(&mut out).unwrap();
            let written = out
                .write_all(b"alignment\n")
                .and_then(|()| display.leave(&mut out));
            assert!(written.is_err());
            // The failure ends the run in the region.
            assert!(cell.alignments.get());
        }
        assert_eq!(sink.routed, routed);
    }
}
