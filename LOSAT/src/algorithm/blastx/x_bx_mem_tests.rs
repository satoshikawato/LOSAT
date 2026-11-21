// No NCBI counterpart (scheduling only): the capped units' task sizes of LOSAT_X_BXCAPFILL.
#[test]
fn capped_units_have_f_tasks_for_each_thread_they_are_left() {
    use super::{filled_subjects_per_task as size, unit_share as share};
    // 50 subjects, K = 8: 16 threads leave each unit 2, so 8 tasks of 7; 64 threads 8, so 32
    // tasks of 2, but a batch of one task 2 again; 6,600 subjects keep the size without the cap.
    assert_eq!(size(4, share(16, 8, false), 50, 64), 7);
    assert_eq!(size(4, share(64, 8, false), 50, 64), 2);
    assert_eq!(size(4, share(64, 8, true), 50, 64), 7);
    assert_eq!(size(4, share(16, 8, false), 6_600, 64), 64);
    assert_eq!(size(4, share(16, 8, false), 6_600, 16), 16);
    // A batch of one task is left the share of 16 threads at K = 8, and never more than others.
    for threads in [9, 16, 17, 32, 48, 64, 128] {
        for k in [2, 4, 8] {
            assert_eq!(share(threads, k, true), share(threads, k, false).min(2));
        }
    }
    // K units then hold a task for every thread they are left, and without the bound of a
    // batch of one task a task for every thread (the rule `unit_cap` keeps without the fill),
    // wherever K units of one subject each would.
    for threads in [9, 16, 17, 32, 48, 64, 128] {
        for one_task in [false, true] {
            for subjects in [2, 10, 50, 200, 543, 6_600] {
                let left = share(threads, 8, one_task);
                let per_task = size(4, left, subjects, 64);
                let tasks = subjects.div_ceil(per_task);
                assert!(
                    tasks >= left || subjects < left,
                    "{threads} threads, {subjects} subjects: {per_task} per task"
                );
                assert!(
                    one_task || 8 * tasks >= threads || 8 * subjects < threads,
                    "{threads} threads, {subjects} subjects: {per_task} per task"
                );
            }
        }
    }
    // One subject at least, also for no subjects.
    assert_eq!(size(4, share(16, 8, false), 0, 64), 1);
    assert_eq!(size(64, share(4_096, 2, false), 1, 16), 1);
}

// No NCBI counterpart (scheduling only). The counts are the run's own (attached to its pool),
// so other tests searching at the same time do not move them.
#[test]
fn a_wave_holds_its_chunks_until_it_ends() {
    let run = super::XCapRun::default();
    let pool = crate::utils::threading::SearchPool::x_serial();
    let pool = run.attach(&pool);
    assert!(!super::in_capped_batch(&pool));
    let batch = super::XCappedBatch::enter(&pool);
    assert!(super::in_capped_batch(&pool));
    // K + 1 less the batch itself, then one while the first wave holds its chunks.
    let first = super::chunk_wave(&pool, 8);
    assert_eq!(first.size(), 8);
    let second = super::chunk_wave(&pool, 8);
    assert_eq!(second.size(), 1);
    drop(first);
    assert_eq!(super::chunk_wave(&pool, 8).size(), 7);
    drop((second, batch));
    assert!(!super::in_capped_batch(&pool));
}
