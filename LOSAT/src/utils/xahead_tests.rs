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
fn ordered_sum(n: usize, window: usize, helpers: usize, leave: bool) -> (u64, Vec<usize>) {
    let compute = |index: usize, scratch: &mut u64| -> u64 {
        *scratch += 1;
        value_of(index)
    };
    let ahead = if leave {
        Ahead::x_new_leaving(window)
    } else {
        Ahead::new(window)
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
    eprintln!("window {window} helpers {helpers} leave {leave}: evaluated ahead {evaluated}, taken {taken}, kept {}", kept.len());
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
            assert_eq!(ordered_sum(5000, window, helpers, false), serial);
            // Helpers that leave once every index is claimed (the owner takes the rest).
            assert_eq!(ordered_sum(5000, window, helpers, true), serial);
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
