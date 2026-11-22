use super::*;
use crate::algorithm::blastp::hsp::BlastpHspList;
use crate::algorithm::blastp::kappa::smith_waterman::SwRedoFailed;

/// A batch status, failed (`failed`, a failed Smith-Waterman redo before) or not.
fn status(failed: bool) -> SwBatchStatus {
    let mut status = SwBatchStatus::default();
    if failed {
        let redone =
            status.redone::<()>(Err(SwRedoFailed.into()), || anyhow::anyhow!("not endless"));
        assert!(matches!(redone, Ok(None)));
    }
    assert_eq!(status.failed(), failed);
    status
}

// Line RG item X4-excl: at a list whose fetch is excluded the Smith-Waterman status is cleared
// unless the early termination test holds (every heap of the batch filled to its cutoff and
// the list's best e-value above 5 times the cutoff); the heaps are read only.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3575
#[test]
fn an_excluded_list_clears_the_status_unless_the_redo_terminates_early() {
    let mark = |best_evalue| super::super::database::FetchExcludedMatch {
        query: 0,
        position: 0,
        best_evalue,
    };
    let empty = [BlastCompoHeap::new(1, 0.002)];
    let mut filled = BlastCompoHeap::new(1, 0.002);
    let list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: Vec::new(),
        best_evalue: 1e-40,
    };
    filled.insert(list, 1e-40, 400, 0);
    assert!(filled.filled_to_cutoff());
    let filled = [filled];
    // A batch of two queries, the second heap empty: the test does not hold.
    let mixed = [filled[0].clone(), BlastCompoHeap::new(1, 0.002)];
    for (heaps, best_evalue, cleared) in [
        (&empty[..], 0.04, true),
        (&filled[..], 0.04, false),
        (&filled[..], 0.01, true),
        (&filled[..], 1e-30, true),
        (&mixed[..], 0.04, true),
    ] {
        for failed_before in [true, false] {
            let mut sw_status = status(failed_before);
            fetch_excluded_reached(&mark(best_evalue), heaps, &mut sw_status);
            let expected = failed_before && !cleared;
            assert_eq!(
                sw_status.failed(),
                expected,
                "{best_evalue} {}",
                heaps.len()
            );
        }
    }
}

#[test]
fn only_excluded_lists_after_the_last_status_count() {
    let mark = |position, best_evalue| super::super::database::FetchExcludedMatch {
        query: 0,
        position,
        best_evalue,
    };
    let empty = [BlastCompoHeap::new(1, 0.002)];
    let mut filled = BlastCompoHeap::new(1, 0.002);
    let list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: Vec::new(),
        best_evalue: 1e-40,
    };
    filled.insert(list, 1e-40, 400, 0);
    let filled = [filled];
    // `position` counts the kept lists before an excluded one: 2 is after the kept list at
    // place 1 and before the one at place 2.
    for (marks, status_position, heaps, cleared) in [
        // Before, at and after the place of the last status: only after it clears.
        (vec![mark(0, 0.04)], Some(1), &empty[..], false),
        (vec![mark(1, 0.04)], Some(1), &empty[..], false),
        (vec![mark(2, 0.04)], Some(1), &empty[..], true),
        (
            vec![mark(0, 0.04), mark(2, 0.04)],
            Some(1),
            &empty[..],
            true,
        ),
        // After it, but the early termination holds: nothing changes.
        (vec![mark(2, 0.04)], Some(1), &filled[..], false),
        (
            vec![mark(2, 0.04), mark(2, 0.005)],
            Some(1),
            &filled[..],
            true,
        ),
        // No redo of the query set the status (a failure of an earlier query of the batch
        // stands): every excluded list counts.
        (vec![mark(0, 0.04)], None, &empty[..], true),
        (vec![mark(0, 0.04)], None, &filled[..], false),
        (vec![mark(0, 0.04), mark(0, 1e-30)], None, &filled[..], true),
        (Vec::new(), None, &empty[..], false),
    ] {
        for failed_before in [true, false] {
            let mut sw_status = status(failed_before);
            query_fetch_excluded_reached(&marks, status_position, heaps, &mut sw_status);
            let expected = failed_before && !cleared;
            assert_eq!(
                sw_status.failed(),
                expected,
                "{status_position:?} {}",
                marks.len()
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3644-3646
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3540
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3560-3575
/// The status of a query batch (`SwBatchStatus`, the one rule of both loops of the
/// Smith-Waterman redo): a redo that ends sets it to 0 and keeps its result, a failed one sets it
/// and keeps nothing, one that does not end gives the caller's error, another error passes as it
/// is; the early termination test leaves it, an excluded fetch after the test sets it to 0; a batch
/// starts at 0.
#[test]
fn the_batch_status_follows_the_last_match_that_the_redo_reaches() {
    let endless = || anyhow::anyhow!("endless here");
    for failed_before in [true, false] {
        let fresh = status(failed_before);
        let mut status = fresh;
        assert_eq!(status.redone(Ok(7), endless).ok(), Some(Some(7)));
        assert!(!status.failed());
        let mut status = fresh;
        let failed = status.redone::<u8>(Err(SwRedoFailed.into()), endless);
        assert_eq!(failed.ok(), Some(None));
        assert!(status.failed());
        let mut status = fresh;
        let error = status
            .redone::<u8>(Err(kappa::smith_waterman::SwRedoEndless.into()), endless)
            .expect_err("a redo that does not end");
        assert_eq!(format!("{error:#}"), "endless here");
        let error = status
            .redone::<u8>(Err(anyhow::anyhow!("a read error")), endless)
            .expect_err("another error");
        assert_eq!(format!("{error:#}"), "a read error");
        // NCBI's test, then its fetch: the test ends a match whether or not the fetch would
        // exclude it, and leaves the status; an excluded fetch clears it; else the redo runs.
        for (terminates_early, excluded, redone, failed_after) in [
            (true, true, false, failed_before),
            (true, false, false, failed_before),
            (false, true, false, false),
            (false, false, true, failed_before),
        ] {
            let mut status = fresh;
            assert_eq!(status.reaches_redo(terminates_early, excluded), redone);
            assert_eq!(status.failed(), failed_after);
        }
    }
    assert!(!SwBatchStatus::default().failed());
}
