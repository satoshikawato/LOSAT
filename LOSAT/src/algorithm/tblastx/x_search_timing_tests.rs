use super::*;

fn stage_sum(times: &XSearchTimes) -> u64 {
    STAGES.iter().map(|(stage, _)| times.stage(*stage)).sum()
}

#[test]
fn the_stage_table_follows_the_discriminants() {
    for (index, (stage, _)) in STAGES.iter().enumerate() {
        assert_eq!(*stage as usize, index, "{stage:?}");
    }
    let mut names: Vec<&str> = STAGES.iter().map(|(_, name)| *name).collect();
    names.sort_unstable();
    names.dedup();
    assert_eq!(names.len(), STAGE_COUNT);
}

#[test]
fn with_the_timer_off_nothing_is_read_or_counted() {
    let times = XSearchTimes::new(false);
    assert!(times.now().is_none());
    {
        let mut lap = times.subject();
        assert!(lap.start.is_none() && lap.last.is_none());
        lap.lap(XStage::Units);
        lap.lap(XStage::Link1);
    }
    assert_eq!(stage_sum(&times), 0);
    assert_eq!(times.subject_ns.load(AtomicOrdering::Relaxed), 0);
    assert_eq!(times.subjects.load(AtomicOrdering::Relaxed), 0);
    assert!(times
        .lines(Instant::now(), None, None, Duration::ZERO)
        .is_none());
}

#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "the Node WASI host can sleep less than asked: a 1 ms sleep timed under 1 ms in 1 run of 6"
)]
fn the_stages_of_a_subject_add_up_to_its_time() {
    let times = XSearchTimes::new(true);
    for _ in 0..3 {
        let mut lap = times.subject();
        std::thread::sleep(Duration::from_millis(2));
        lap.lap(XStage::Translate);
        std::thread::sleep(Duration::from_millis(3));
        lap.lap(XStage::Link1);
        std::thread::sleep(Duration::from_millis(1));
    }
    let total = times.subject_ns.load(AtomicOrdering::Relaxed);
    assert_eq!(stage_sum(&times), total);
    assert_eq!(times.subjects.load(AtomicOrdering::Relaxed), 3);
    assert!(times.stage(XStage::Translate) >= 6_000_000);
    assert!(times.stage(XStage::Link1) >= 9_000_000);
    assert!(times.stage(XStage::End) >= 3_000_000);
    assert_eq!(times.stage(XStage::Units), 0);
}

#[test]
fn a_subject_that_returns_early_is_all_in_the_last_stage() {
    let times = XSearchTimes::new(true);
    drop(times.subject());
    assert_eq!(times.subjects.load(AtomicOrdering::Relaxed), 1);
    assert_eq!(
        times.stage(XStage::End),
        times.subject_ns.load(AtomicOrdering::Relaxed)
    );
}

#[test]
fn the_first_line_adds_up_to_search_total() {
    let times = XSearchTimes::new(true);
    let search_start = Instant::now();
    let subjects_start = search_start + Duration::from_millis(250);
    let subjects_end = subjects_start + Duration::from_millis(1500);
    let [search, subject] = times
        .lines(
            search_start,
            Some(subjects_start),
            Some(subjects_end),
            Duration::from_millis(2000),
        )
        .expect("timer on");
    assert_eq!(
        search,
        "[TIMING] x_search: setup 0.250s subjects 1.500s collect 0.250s"
    );
    assert!(subject.starts_with("[TIMING] x_subject: enter 0.000s translate 0.000s "));
    assert!(subject.ends_with(" end 0.000s sum 0.000s (subjects=0)"));
    assert_eq!(subject.matches(" 0.000s").count(), STAGE_COUNT + 1);
}

#[test]
fn the_counters_are_shared_between_the_threads_of_a_search() {
    fn shared<S: Send + Sync>() {}
    shared::<XSearchTimes>();
}
