// No NCBI counterpart: NCBI frees its tables after the subject loop and has no unit test of
// an earlier release. These tests check the ownership rules of the module doc on a payload
// that counts its drops.
use super::*;

/// A stand-in for the tables: counts how often it is dropped.
struct Payload {
    drops: Arc<AtomicUsize>,
    value: u32,
}

impl Drop for Payload {
    fn drop(&mut self) {
        self.drops.fetch_add(1, AtomicOrdering::SeqCst);
    }
}

fn new_state(stages: usize) -> (XSeedState<Payload>, Arc<AtomicUsize>) {
    let drops = Arc::new(AtomicUsize::new(0));
    let payload = Payload {
        drops: Arc::clone(&drops),
        value: 7,
    };
    (XSeedState::new(payload, stages, None), drops)
}

fn freed(drops: &AtomicUsize) -> usize {
    drops.load(AtomicOrdering::SeqCst)
}

#[test]
fn one_subject_holds_the_tables_until_its_stage_ends() {
    let (state, drops) = new_state(1);
    let stage = state.enter(0);
    assert_eq!(stage.tables().value, 7);
    assert_eq!(freed(&drops), 0);
    assert!(stage.end());
    assert_eq!(freed(&drops), 1);
    drop(state);
    assert_eq!(freed(&drops), 1);
}

#[test]
fn subjects_in_order_free_the_tables_at_the_last_end() {
    let (state, drops) = new_state(4);
    for s_idx in 0..4 {
        let stage = state.enter(s_idx);
        assert_eq!(stage.tables().value, 7);
        assert_eq!(freed(&drops), 0);
        assert_eq!(stage.end(), s_idx == 3, "subject {s_idx}");
        assert_eq!(freed(&drops), usize::from(s_idx == 3), "subject {s_idx}");
    }
}

#[test]
fn a_stage_dropped_without_end_counts() {
    // Not the last one: the next stage's `end` is the release.
    let (state, drops) = new_state(2);
    drop(state.enter(0));
    assert_eq!(freed(&drops), 0);
    assert!(state.enter(1).end());
    assert_eq!(freed(&drops), 1);

    // The last one: its drop is the release.
    let (state, drops) = new_state(2);
    assert!(!state.enter(0).end());
    assert_eq!(freed(&drops), 0);
    {
        let stage = state.enter(1);
        assert_eq!(stage.tables().value, 7);
        assert_eq!(freed(&drops), 0);
    }
    assert_eq!(freed(&drops), 1);
}

#[test]
#[cfg_attr(panic = "abort", ignore = "catches a panic: the build aborts on panic")]
fn a_stage_dropped_by_a_panic_counts() {
    let (state, drops) = new_state(1);
    let unwound = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let stage = state.enter(0);
        assert_eq!(stage.tables().value, 7);
        panic!("a panic inside the seed stage");
    }));
    assert!(unwound.is_err());
    assert_eq!(freed(&drops), 1);
}

#[test]
fn stages_open_at_once_free_the_tables_once_after_the_last() {
    let (state, drops) = new_state(3);
    let a = state.enter(0);
    let b = state.enter(1);
    let c = state.enter(2);
    assert!(!c.end());
    assert_eq!(freed(&drops), 0);
    assert!(!a.end());
    assert_eq!(freed(&drops), 0);
    assert_eq!(b.tables().value, 7);
    assert!(b.end());
    assert_eq!(freed(&drops), 1);
    drop(state);
    assert_eq!(freed(&drops), 1);
}

#[cfg(any(not(target_arch = "wasm32"), feature = "wasm-threads"))]
#[test]
fn stages_on_several_threads_free_the_tables_once_after_the_last() {
    const THREADS: usize = 8;
    for round in 0..200 {
        let (state, drops) = new_state(THREADS);
        // Even rounds: every stage is open before the first one ends. Odd rounds: the
        // threads enter and end as they come.
        let all_open = std::sync::Barrier::new(if round % 2 == 0 { THREADS } else { 1 });
        let released = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            for s_idx in 0..THREADS {
                let (state, drops, all_open, released) = (&state, &drops, &all_open, &released);
                scope.spawn(move || {
                    let stage = state.enter(s_idx);
                    all_open.wait();
                    assert_eq!(stage.tables().value, 7);
                    // This stage is open, so the tables are alive.
                    assert_eq!(freed(drops), 0);
                    if stage.end() {
                        // The release is the last end: every other stage has ended.
                        assert_eq!(freed(drops), 1);
                        released.fetch_add(1, AtomicOrdering::SeqCst);
                    }
                });
            }
        });
        assert_eq!(released.load(AtomicOrdering::SeqCst), 1, "round {round}");
        assert_eq!(freed(&drops), 1, "round {round}");
    }
}

#[test]
#[should_panic(
    expected = "the seed-stage tables of the batch were released before the seed stage of subject 5 started"
)]
fn enter_after_the_release_stops_the_run() {
    let (state, _drops) = new_state(1);
    assert!(state.enter(0).end());
    let _stage = state.enter(5);
}

#[test]
fn a_count_that_is_too_low_keeps_the_tables_until_the_last_open_stage_ends() {
    // One stage announced, two open at once.
    let (state, drops) = new_state(1);
    let a = state.enter(0);
    let b = state.enter(1);
    // The announced count is used up: the batch's hold goes, the tables stay.
    assert!(a.end());
    assert_eq!(freed(&drops), 0);
    assert_eq!(b.tables().value, 7);
    // The count stays at zero; the last open stage frees the tables and reports no release.
    assert!(!b.end());
    assert_eq!(freed(&drops), 1);
}

#[test]
fn a_count_that_is_too_high_keeps_the_tables_until_the_state_is_dropped() {
    let (state, drops) = new_state(3);
    assert!(!state.enter(0).end());
    assert!(!state.enter(1).end());
    assert_eq!(freed(&drops), 0);
    drop(state);
    assert_eq!(freed(&drops), 1);
}

#[test]
fn a_batch_without_subjects_keeps_the_tables_until_the_state_is_dropped() {
    let (state, drops) = new_state(0);
    assert_eq!(freed(&drops), 0);
    drop(state);
    assert_eq!(freed(&drops), 1);
}

#[test]
fn workers_get_no_diagonal_table_only_when_the_batch_is_units_only() {
    let cells = 1usize << 21;
    for (units_only, expected) in [(true, 0), (false, cells)] {
        assert_eq!(
            worker_table_cells(units_only, cells),
            expected,
            "units_only {units_only}"
        );
    }
}

#[test]
fn drop_worker_table_leaves_no_allocation() {
    let mut table = vec![0u32; 1024];
    drop_worker_table(&mut table);
    assert!(table.is_empty());
    assert_eq!(table.capacity(), 0);
    // A table of zero cells (no reference frame loop in the batch) stays as it is.
    drop_worker_table(&mut table);
    assert_eq!(table.capacity(), 0);
}

#[test]
fn the_release_adds_one_to_the_counter_of_its_state() {
    fn counter() -> &'static AtomicU64 {
        Box::leak(Box::new(AtomicU64::new(0)))
    }
    let count = |c: &AtomicU64| c.load(AtomicOrdering::SeqCst);

    // Once, at the release, also when further stages were open.
    let released = counter();
    let state = XSeedState::new(7u32, 2, Some(released));
    let a = state.enter(0);
    let b = state.enter(1);
    let c = state.enter(2);
    assert!(!a.end());
    assert_eq!(count(released), 0);
    assert!(b.end());
    assert_eq!(count(released), 1);
    assert!(!c.end());
    drop(state);
    assert_eq!(count(released), 1);
}

/// A lookup table without entries (every field of the type is public).
fn empty_lookup() -> BlastAaLookupTable {
    BlastAaLookupTable {
        backbone: Vec::new(),
        overflow: Vec::new(),
        pv: Vec::new(),
        frame_bases: Vec::new(),
        num_contexts: 0,
        query_length: 0,
        word_length: 3,
        alphabet_size: 28,
        charsize: 5,
        mask: 0,
        longest_chain: 11,
        threshold: 13,
        row_max: Vec::new(),
    }
}

#[test]
fn the_two_holders_of_a_batch_count_their_stages_apart() {
    // A batch without a unit path (the reference frame loop takes every subject): its
    // workers keep their diagonal table.
    let batch = XSeedBatch::<Blosum62Score>::new(empty_lookup(), None, 2);
    assert_eq!(batch.worker_table_cells(64), 64);

    // First subject: the units stage ends first, then the lookup stage.
    let lookup_stage = batch.lookup.enter(0);
    let units_stage = batch.units.enter(0);
    assert!(units_stage.tables().is_none());
    assert!(!units_stage.end());
    assert_eq!(lookup_stage.tables().longest_chain, 11);
    assert!(!lookup_stage.end());

    // Last subject: the end of its units stage releases the unit holder only; the lookup
    // table is still read through the stage that is open.
    let lookup_stage = batch.lookup.enter(1);
    let units_stage = batch.units.enter(1);
    assert!(units_stage.end());
    assert_eq!(lookup_stage.tables().threshold, 13);
    assert!(lookup_stage.end());
}

#[test]
#[should_panic(
    expected = "the seed-stage tables of the batch were released before the seed stage of subject 1 started"
)]
fn a_unit_stage_after_the_release_of_the_unit_holder_stops_the_run() {
    let batch = XSeedBatch::<Blosum62Score>::new(empty_lookup(), None, 1);
    let lookup_stage = batch.lookup.enter(0);
    assert!(batch.units.enter(0).end());
    // The lookup holder still serves the open stage, and would serve a new one.
    assert_eq!(lookup_stage.tables().longest_chain, 11);
    let _units_stage = batch.units.enter(1);
}

#[test]
fn the_holders_are_shared_between_the_threads_of_a_search() {
    fn shared<S: Send + Sync>() {}
    shared::<XSeedBatch>();
}
