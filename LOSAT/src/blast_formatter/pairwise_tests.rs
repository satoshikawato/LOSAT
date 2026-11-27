use super::*;
use crate::report::archive::read::options::SeqRange;

fn options(program: EProgram, setters: Vec<Setter>) -> BuiltOptions {
    BuiltOptions {
        program,
        setters,
        entrez_query: None,
        first_db_seq: None,
        final_db_seq: None,
        gi_list: None,
        negative_gi_list: None,
        taxid_list: None,
        negative_taxid_list: None,
        db_filtering_algorithm_id: None,
        db_filtering_algorithm_key: None,
        subject_masking_type: None,
        query_masks: None,
        query_range: SeqRange::EMPTY,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:136-219
// The factory defaults of each program's handle, which the epilog prints where the archive
// sets nothing (set A: blastn 2 -3 5 2, megablast 1 -2 0 0, dc-megablast window 40, rpsblast
// without a word threshold, the protein handles' 11 1 11 40).
#[test]
fn epilog_values_are_the_handles_defaults() {
    let values = |program| {
        let v = EpilogValues::of(&options(program, Vec::new()));
        (
            v.reward,
            v.penalty,
            v.gap_open,
            v.gap_extend,
            v.word_threshold,
            v.window_size,
        )
    };
    assert_eq!(values(EProgram::Blastn), (2, -3, 5, 2, 0, 0));
    assert_eq!(values(EProgram::Megablast), (1, -2, 0, 0, 0, 0));
    assert_eq!(values(EProgram::DiscMegablast), (2, -3, 5, 2, 0, 40));
    let (_, _, open, extend, threshold, window) = values(EProgram::RpsBlast);
    assert_eq!((open, extend, threshold, window), (11, 1, 0, 40));
    for program in [
        EProgram::Blastp,
        EProgram::PsiBlast,
        EProgram::DeltaBlast,
        EProgram::PhiBlastp,
    ] {
        let (_, _, open, extend, threshold, window) = values(program);
        assert_eq!(
            (open, extend, threshold, window),
            (11, 1, 11, 40),
            "{program:?}"
        );
    }
}

// The archive's calls on the factory options, in order, the last one of each kind kept
// (blastp-short's PAM30 9 1 16 15 of set A; an ungapped search; sum statistics).
#[test]
fn epilog_values_take_the_archives_calls() {
    let v = EpilogValues::of(&options(
        EProgram::Blastp,
        vec![
            Setter::GapOpeningCost(10),
            Setter::GapOpeningCost(9),
            Setter::WordThreshold(16),
            Setter::WindowSize(15),
            Setter::GappedMode(false),
            Setter::SumStatisticsMode(true),
        ],
    ));
    assert_eq!(
        (
            v.gap_open,
            v.gap_extend,
            v.word_threshold,
            v.window_size,
            v.gapped,
            v.sum_statistics
        ),
        (9, 1, 16, 15, false, true)
    );
    let v = EpilogValues::of(&options(
        EProgram::Blastn,
        vec![Setter::MatchReward(1), Setter::MismatchPenalty(-3)],
    ));
    assert_eq!((v.reward, v.penalty, v.gapped), (1, -3, true));
}

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
