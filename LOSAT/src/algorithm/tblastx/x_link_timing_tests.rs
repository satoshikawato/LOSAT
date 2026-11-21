use super::*;
use std::time::Duration;

fn load(counter: &AtomicU64) -> u64 {
    counter.load(AtomicOrdering::Relaxed)
}

fn counters(timer: &XLinkTimer, call: LinkCall) -> &CallTimes {
    &timer.times.as_ref().expect("timer on").calls[call as usize]
}

/// One timed call that spends a little time in every stage and in two per-group calls.
fn one_call(timer: &XLinkTimer, call: LinkCall) {
    timer.arm(call);
    {
        let mut lap = LinkLap::start();
        for (stage, _) in STAGES {
            if stage == LinkStage::Kernels {
                lap.sizes(10, [4usize, 6].into_iter());
                let kernels = lap.kernels();
                for millis in [1u64, 3] {
                    kernels.time(|| std::thread::sleep(Duration::from_millis(millis)));
                }
                lap.add_kernels(&kernels);
            }
            std::thread::sleep(Duration::from_millis(1));
            lap.lap(stage);
        }
    }
    timer.disarm();
}

// No NCBI counterpart (a timer).
#[test]
fn the_stage_table_follows_the_discriminants() {
    for (index, (stage, _)) in STAGES.iter().enumerate() {
        assert_eq!(*stage as usize, index);
    }
    for (index, (call, _)) in CALLS.iter().enumerate() {
        assert_eq!(*call as usize, index);
    }
}

// No NCBI counterpart (a timer). Over several calls of a batch the sizes and per-group sums
// add up, while `largest_group` and `kernel_max` keep the largest of the calls.
#[test]
fn the_largest_group_and_kernel_are_the_largest_over_the_calls() {
    let timer = XLinkTimer::new(true);
    for (sizes, kernel_max) in [([4usize, 6], 7u64), ([2usize, 3], 9), ([5usize, 1], 2)] {
        timer.arm(LinkCall::First);
        {
            let lap = LinkLap::start();
            lap.sizes(sizes.iter().sum(), sizes.into_iter());
            let kernels = KernelTimes {
                enabled: true,
                sum_ns: AtomicU64::new(kernel_max + 1),
                max_ns: AtomicU64::new(kernel_max),
            };
            lap.add_kernels(&kernels);
        }
        timer.disarm();
    }
    let first = counters(&timer, LinkCall::First);
    assert_eq!(
        (
            load(&first.hsps),
            load(&first.groups),
            load(&first.largest_group)
        ),
        (21, 6, 6)
    );
    assert_eq!(
        (load(&first.kernel_sum_ns), load(&first.kernel_max_ns)),
        (21, 9)
    );
}

// No NCBI counterpart (a timer). Off: `arm` sets no flag, so the call's clock measures
// nothing, its per-group timer reads no clock, and there are no lines.
#[test]
fn with_the_timer_off_nothing_is_armed_or_counted() {
    let timer = XLinkTimer::new(false);
    timer.arm(LinkCall::First);
    let mut lap = LinkLap::start();
    assert!(lap.running.is_none());
    lap.lap(LinkStage::Sort);
    lap.sizes(3, [3usize].into_iter());
    let kernels = lap.kernels();
    assert_eq!(kernels.time(|| 7), 7);
    assert_eq!((load(&kernels.sum_ns), load(&kernels.max_ns)), (0, 0));
    lap.add_kernels(&kernels);
    drop(lap);
    timer.disarm();
    assert!(timer.lines().is_none());
    assert!(LinkLap::off().running.is_none());
}

// No NCBI counterpart (a timer). The stages of a call add up to its whole time, the call is
// counted under the kind it was armed with, and the sizes and per-group times are its own.
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "the Node WASI host can sleep less than asked: a 1 ms sleep timed under 1 ms in 1 run of 6"
)]
fn the_stages_of_a_call_add_up_to_its_time() {
    let timer = XLinkTimer::new(true);
    one_call(&timer, LinkCall::Second);
    let second = counters(&timer, LinkCall::Second);
    assert_eq!(load(&second.calls), 1);
    let stage_sum: u64 = second.stage_ns.iter().map(load).sum();
    let call = load(&second.call_ns);
    assert!(second.stage_ns.iter().all(|stage| load(stage) >= 1_000_000));
    assert_eq!(stage_sum, call);
    assert_eq!(
        (
            load(&second.hsps),
            load(&second.groups),
            load(&second.largest_group)
        ),
        (10, 2, 6)
    );
    let (sum, max) = (load(&second.kernel_sum_ns), load(&second.kernel_max_ns));
    assert!(max >= 3_000_000 && sum >= max + 1_000_000, "{sum} {max}");
    assert!(sum <= load(&second.stage_ns[LinkStage::Kernels as usize]));

    let first = counters(&timer, LinkCall::First);
    assert_eq!(load(&first.calls), 0);
    assert_eq!(load(&first.call_ns), 0);
    let lines = timer.lines().expect("timer on");
    assert!(lines[0].starts_with("[TIMING] x_link1: reset 0.000s sort 0.000s groups 0.000s kernels 0.000s assemble 0.000s final_sorts 0.000s replay 0.000s sum 0.000s (calls=0 hsps=0 groups=0 largest_group=0 kernel_sum 0.000s kernel_max 0.000s)"), "{}", lines[0]);
    // The sleeps above have no upper bound (several milliseconds each on a loaded machine),
    // so the times of this line are not compared with text.
    assert!(
        lines[1].starts_with("[TIMING] x_link2: reset "),
        "{}",
        lines[1]
    );
    assert!(
        lines[1].contains("s sum ")
            && lines[1].contains("s (calls=1 hsps=10 groups=2 largest_group=6 kernel_sum ")
            && lines[1].contains("s kernel_max ")
            && lines[1].ends_with("s)"),
        "{}",
        lines[1]
    );
}

// No NCBI counterpart (a timer). The flag is for one call: a second call of the thread
// without `arm` is not timed, and a flag that no call took is cleared by `disarm`.
#[test]
fn the_flag_is_taken_once_and_cleared_when_not_taken() {
    let timer = XLinkTimer::new(true);
    timer.arm(LinkCall::First);
    let first = LinkLap::start();
    assert!(first.running.is_some());
    let untimed = LinkLap::start();
    assert!(untimed.running.is_none());
    drop(untimed);
    drop(first);
    assert_eq!(load(&counters(&timer, LinkCall::First).calls), 1);

    timer.arm(LinkCall::Second);
    timer.disarm();
    assert!(LinkLap::start().running.is_none());
    assert_eq!(load(&counters(&timer, LinkCall::Second).calls), 0);
}

// No NCBI counterpart (a timer). A thread that is given another subject while its own
// linking call waits for the pool: the inner calls arm, take and clear the flag before the
// outer call ends, and every call is counted once, under its own batch and kind.
#[test]
fn a_call_nested_in_a_running_call_is_counted_separately() {
    let outer_timer = XLinkTimer::new(true);
    let inner_timer = XLinkTimer::new(true);
    outer_timer.arm(LinkCall::First);
    {
        let mut outer = LinkLap::start();
        outer.lap(LinkStage::Groups);
        one_call(&inner_timer, LinkCall::First);
        one_call(&inner_timer, LinkCall::Second);
        outer.lap(LinkStage::Kernels);
    }
    outer_timer.disarm();
    assert_eq!(load(&counters(&outer_timer, LinkCall::First).calls), 1);
    assert_eq!(load(&counters(&outer_timer, LinkCall::Second).calls), 0);
    assert_eq!(load(&counters(&inner_timer, LinkCall::First).calls), 1);
    assert_eq!(load(&counters(&inner_timer, LinkCall::Second).calls), 1);
    let outer_kernels =
        load(&counters(&outer_timer, LinkCall::First).stage_ns[LinkStage::Kernels as usize]);
    let inner: u64 = [LinkCall::First, LinkCall::Second]
        .iter()
        .map(|call| load(&counters(&inner_timer, *call).call_ns))
        .sum();
    assert!(outer_kernels >= inner, "{outer_kernels} {inner}");
}

// No NCBI counterpart (a timer). The counters are shared by the threads of a search, and the
// per-group timer by the threads of the pool.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn the_counters_are_shared_between_threads() {
    fn shared<S: Send + Sync>() {}
    shared::<XLinkTimer>();
    shared::<KernelTimes>();
    let timer = XLinkTimer::new(true);
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| one_call(&timer, LinkCall::First));
        }
    });
    let first = counters(&timer, LinkCall::First);
    assert_eq!(load(&first.calls), 4);
    assert_eq!(load(&first.hsps), 40);
    let stage_sum: u64 = first.stage_ns.iter().map(load).sum();
    assert_eq!(stage_sum, load(&first.call_ns));
}
