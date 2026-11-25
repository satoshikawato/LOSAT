use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919
// The tests below run a stand-in for this ordered loop: a value that depends on the
// index alone, consumed in order.
// The owner's result must not depend on what the helpers do.
#[cfg(feature = "parallel")]
fn value_of(index: usize) -> u64 {
    let mut x = index as u64 ^ 0x9E37_79B9_7F4A_7C15;
    for _ in 0..(index % 7) * 3000 {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
    }
    x
}

#[cfg(feature = "parallel")]
fn ordered_sum(n: usize, window: usize, helpers: usize, form: &str) -> (u64, Vec<usize>) {
    let compute = |index: usize, scratch: &mut u64| -> u64 {
        *scratch += 1;
        value_of(index)
    };
    let ahead = match form {
        "leaving" => Ahead::x_new_leaving(window),
        "napping" => Ahead::x_new_napping(window),
        _ => Ahead::new(window),
    };
    let mut kept: Vec<usize> = Vec::new();
    let mut sum = 0u64;
    let mut owner_scratch = 0u64;
    let scratches: Vec<u64> = vec![0; helpers];
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(helpers + 1)
        .build()
        .unwrap();
    pool.install(|| {
        with_helpers(Some(&ahead), scratches, &compute, || {
            for index in 0..n {
                // What the owner knows changes with every value it keeps.
                let live = |j: usize, sum: u64| !(j as u64).wrapping_add(sum).is_multiple_of(3);
                ahead.announce(index, n, |j| live(j, sum));
                if !live(index, sum) {
                    ahead.skip(index);
                    continue;
                }
                let value = ahead
                    .take(index, &mut owner_scratch, &compute)
                    .unwrap_or_else(|| compute(index, &mut owner_scratch));
                sum = sum.wrapping_mul(31).wrapping_add(value);
                kept.push(index);
            }
        })
    });
    let (evaluated, taken) = ahead.counts();
    assert!(taken <= evaluated);
    eprintln!("window {window} helpers {helpers} form {form}: evaluated ahead {evaluated}, taken {taken}, kept {}", kept.len());
    (sum, kept)
}

#[cfg(feature = "parallel")]
#[test]
fn x_ahead_owner_result_is_independent_of_helpers() {
    let serial = {
        let mut kept = Vec::new();
        let mut sum = 0u64;
        for index in 0..5000usize {
            if (index as u64).wrapping_add(sum).is_multiple_of(3) {
                continue;
            }
            sum = sum.wrapping_mul(31).wrapping_add(value_of(index));
            kept.push(index);
        }
        (sum, kept)
    };
    for &(window, helpers) in &[(2, 1), (8, 1), (64, 1), (64, 3), (5, 2), (256, 7)] {
        for _ in 0..3 {
            assert_eq!(ordered_sum(5000, window, helpers, "staying"), serial);
            // Helpers that leave once every index is claimed (the owner takes the rest).
            assert_eq!(ordered_sum(5000, window, helpers, "leaving"), serial);
            // And whose waits nap (BLASTX's composition redo).
            assert_eq!(ordered_sum(5000, window, helpers, "napping"), serial);
        }
    }
}

#[test]
fn x_ahead_without_pool_is_the_plain_loop() {
    let compute = |index: usize, _: &mut ()| index * 2;
    let ahead = Ahead::new(8);
    let total = with_helpers(Some(&ahead), vec![(), ()], &compute, || {
        let mut total = 0;
        for index in 0..100 {
            ahead.announce(index, 100, |_| true);
            total += ahead
                .take(index, &mut (), &compute)
                .unwrap_or_else(|| compute(index, &mut ()));
        }
        total
    });
    assert_eq!(total, 9900);
    assert_eq!(ahead.counts(), (0, 0));
}

// No NCBI counterpart (scheduling only): a guess of "not yet" (`x_announce_while`, BLASTX's
// composition redo) stops the offers at its index; a later call starts there, and the helpers that
// leave once every index is claimed stay while the offers are stopped.
#[test]
fn x_ahead_stopped_offers_resume_where_they_stopped() {
    let ahead: Ahead<usize> = Ahead::new(8);
    let mut asked = Vec::new();
    ahead.x_announce_while(0, 20, |j| {
        asked.push(j);
        (j < 3).then_some(true)
    });
    assert_eq!(asked, [0, 1, 2, 3]);
    assert_eq!(ahead.announced.0.load(Ordering::Relaxed), 3);
    asked.clear();
    // The owner at 1: the window [1, 9) from the stopped index 3 on.
    ahead.x_announce_while(1, 20, |j| {
        asked.push(j);
        Some(j != 5)
    });
    assert_eq!(asked, [3, 4, 5, 6, 7, 8]);
    assert_eq!(ahead.announced.0.load(Ordering::Relaxed), 9);
    let claimed: Vec<usize> = std::iter::from_fn(|| ahead.claim()).collect();
    assert_eq!(claimed, (0..9).collect::<Vec<_>>());
    // Index 5 was closed: nobody evaluates it; the others are on offer.
    let compute = |index: usize, _: &mut ()| index * 10;
    ahead.evaluate(5, &mut (), &compute);
    ahead.evaluate(4, &mut (), &compute);
    assert_eq!(ahead.take(5, &mut (), &compute), None);
    assert_eq!(ahead.take(4, &mut (), &compute), Some(40));
    // `announce` and `x_announce_until` are the same loop: `until` stops the offers there.
    let ahead: Ahead<usize> = Ahead::new(8);
    ahead.x_announce_until(0, 2, 20, |_| true);
    assert_eq!(ahead.announced.0.load(Ordering::Relaxed), 2);
    ahead.announce(0, 20, |_| true);
    assert_eq!(ahead.announced.0.load(Ordering::Relaxed), 8);
}

#[cfg(feature = "parallel")]
#[test]
fn x_ahead_helpers_stay_while_the_offers_are_stopped() {
    let wait_until = |done: &dyn Fn() -> bool, what: &str| {
        let start = std::time::Instant::now();
        while !done() {
            assert!(
                start.elapsed() < std::time::Duration::from_secs(60),
                "timed out waiting for {what}"
            );
            std::thread::yield_now();
        }
    };
    let compute = |index: usize, _: &mut ()| index + 100;
    // Helpers that leave once every index is claimed, and nap (BLASTX's composition redo).
    let ahead: Ahead<usize> = Ahead::x_new_napping(8);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    pool.install(|| {
        with_helpers(Some(&ahead), [()], &compute, || {
            // Offers 0 and 1, then "not yet" at 2 of a loop of 6.
            ahead.x_announce_while(0, 6, |j| (j < 2).then_some(true));
            wait_until(&|| ahead.counts().0 == 2, "the helper's first two values");
            // Every offered index is claimed; the helper must not leave.
            std::thread::sleep(std::time::Duration::from_millis(20));
            ahead.x_announce_while(2, 6, |_| Some(true));
            // The owner claims nothing here: only the helper evaluates the later four.
            wait_until(&|| ahead.counts().0 == 6, "the helper's later values");
            for index in 0..6 {
                assert_eq!(ahead.take(index, &mut (), &compute), Some(index + 100));
            }
        })
    });
    assert_eq!(ahead.counts(), (6, 6));
}
