use super::*;

#[test]
fn the_line_names_every_phase_in_order() {
    let l = line(
        Some(Duration::from_millis(1500)),
        Duration::from_millis(1000),
    );
    assert!(l.starts_with("[TIMING] common: since_start 1.500s outside_search 0.500s read_fasta "));
    assert!(l.contains(" search 1.000s output "), "{l}");
    let names: Vec<&str> = l.split_whitespace().skip(2).step_by(2).collect::<Vec<_>>();
    assert_eq!(
        names,
        [
            "since_start",
            "outside_search",
            "read_fasta",
            "encode",
            "translate",
            "seg",
            "build_lookup",
            "search",
            "output",
            "(search_sum"
        ]
    );
}

// Batches searched side by side count once in `search`, and a gap between intervals is
// outside it.
#[test]
fn the_search_is_the_union_of_its_intervals() {
    let t = Instant::now();
    let at = |ms: u64| t + Duration::from_millis(ms);
    let mut spans = vec![
        (at(10), at(40)),
        (at(0), at(20)),
        (at(30), at(35)),
        (at(50), at(60)),
        (at(60), at(70)),
    ];
    assert_eq!(union(&mut spans), Duration::from_millis(60));
    assert_eq!(union(&mut []), Duration::ZERO);
}

#[test]
fn durations_add_up_per_phase() {
    let before = NS[Phase::Seg as usize].load(Ordering::Relaxed);
    add_duration(Phase::Seg, Duration::from_nanos(7));
    add_duration(Phase::Seg, Duration::from_nanos(5));
    assert!(NS[Phase::Seg as usize].load(Ordering::Relaxed) >= before + 12);
    add(Phase::Seg, None);
}
