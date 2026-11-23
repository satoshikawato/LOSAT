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

// No NCBI counterpart (scheduling only): a reader batch searched whole is left the threads of its
// places in the batch flow without the cap (one a place, at most the pool's), its subjects as f
// tasks for each of them; a unit of one place, and a search without a pool, keep their schedule.
#[test]
fn a_batch_searched_whole_is_left_the_threads_of_its_places() {
    let options = crate::algorithm::blastx::stage_c_tests::options();
    let serial = crate::utils::threading::SearchPool::x_serial();
    assert_eq!(super::whole_batch_threads(&serial, &options, 40, 10), None);
    crate::utils::threading::with_search_pool(4, "blastx", |pool| {
        // The run's own counts: no capped batch of another test is searched in it.
        let run = super::XCapRun::default();
        let pool = run.attach(pool);
        let share = |places| super::whole_batch_threads(&pool, &options, 40, places);
        assert_eq!((share(0), share(1)), (None, None));
        assert_eq!((share(2), share(3), share(4)), (Some(2), Some(3), Some(4)));
        assert_eq!((share(10), share(usize::MAX)), (Some(4), Some(4)));
        for threads in [2, 4] {
            assert_eq!(
                super::flow_threads(&pool, &options, 40, Some(threads)),
                threads
            );
            let per_task = super::flow_subjects_per_task(&pool, &options, 40, 64, Some(threads));
            assert_eq!(per_task, super::capped_subjects_per_task(threads, 40, 64));
        }
        // With the built-in f = 4: 40 subjects on 4 threads in tasks of 3 (14 tasks), on 2 in
        // tasks of 5; a database of many subjects keeps tasks of 64.
        if super::cap_fill() == Some(4) {
            let per_task = |subjects, threads| {
                super::flow_subjects_per_task(&pool, &options, subjects, 64, Some(threads))
            };
            assert_eq!((per_task(40, 4), per_task(40, 2)), (3, 5));
            assert_eq!(per_task(100_000, 4), 64);
            // At most ceil(subjects / 4) threads: 8 subjects in tasks of one on 2 threads, 1 or
            // none on 1.
            let threads = |subjects| super::flow_threads(&pool, &options, subjects, Some(4));
            assert_eq!((threads(8), threads(9), threads(16)), (2, 3, 4));
            assert_eq!((threads(1), threads(0)), (1, 1));
            assert_eq!(per_task(8, 4), 1);
        }
        // A unit of one place: the pool's threads and tasks of 64, as before.
        assert_eq!(super::flow_threads(&pool, &options, 40, None), 4);
        assert_eq!(
            super::flow_subjects_per_task(&pool, &options, 40, 64, None),
            64
        );
        Ok(())
    })
    .unwrap();
}
