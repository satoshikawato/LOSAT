use super::*;

// ---- EXPERIMENT (LOSAT_X_BPDIAG): the batched cell test against the per-hit test ----

/// One call of the candidate closure: query offset, subject offset, last_hit, diff, and the
/// entry before the call.
type Call = (i32, i32, i32, i32, u32);

fn next(seed: &mut u64) -> u32 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 32) as u32
}

fn entry_of(raw: u32) -> DiagStruct {
    let mut entry = DiagStruct::default();
    entry.set_raw_bits(raw);
    entry
}

/// The modes this machine runs: the per-hit test (the reference), the scalar batch and, with
/// AVX2, the AVX2 batches (blastp's and blastx's, LOSAT_X_BXCROSS's among them where the CPU
/// also has POPCNT, with its calls deferred, LOSAT_X_BXCROSSDEFER, and not).
fn cell_modes() -> Vec<XCellTest> {
    #[allow(unused_mut)]
    let mut modes = vec![XDiagMode::PerHit, XDiagMode::Scalar];
    #[cfg(target_arch = "x86_64")]
    if x_avx2_available() {
        modes.push(XDiagMode::Avx2);
        modes.push(XDiagMode::Avx2Defer);
        modes.push(XDiagMode::Avx2DeferMasked);
    } else {
        eprintln!("no AVX2 on this CPU: the AVX2 batch is not compared");
    }
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    if x_cross128_available() {
        modes.push(XDiagMode::W128);
    }
    #[allow(unused_mut)]
    let mut modes: Vec<XCellTest> = modes.into_iter().map(XCellTest::Diag).collect();
    #[cfg(target_arch = "x86_64")]
    if x_cross_available() {
        modes.push(XCellTest::Cross { defer: false });
        modes.push(XCellTest::Cross { defer: true });
    }
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    if x_cross128_available() {
        modes.push(XCellTest::Cross128);
    }
    modes
}

/// The per-hit test, the reference of the other cell tests.
const PER_HIT: XCellTest = XCellTest::Diag(XDiagMode::PerHit);

/// The runs of [`cell_modes`] with the LOSAT_X_BXCROSS chunk each takes: LOSAT_X_BXCROSS's
/// modes with the chunk of the test hook and with the built-in one, the other modes once (they
/// take no chunk; 0).
fn cell_runs() -> Vec<(XCellTest, usize)> {
    let mut runs = Vec::new();
    for mode in cell_modes() {
        #[cfg(target_arch = "x86_64")]
        if let XCellTest::Cross { .. } = mode {
            runs.push((mode, X_CROSS_TEST_CHUNK));
            runs.push((mode, X_CROSS_CHUNK));
            continue;
        }
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        if mode == XCellTest::Cross128 {
            runs.push((mode, X_CROSS_TEST_CHUNK));
            runs.push((mode, X_CROSS_CHUNK));
            continue;
        }
        runs.push((mode, 0));
    }
    runs
}

/// The candidate closure of the tests: writes the candidate's entry as the extension does
/// (with the flag for some), so a batch that stored its entries after the call would leave
/// another table.
fn extend(query_offset: i32, subject_offset: i32, entry: &mut DiagStruct, diag_offset: i32) {
    if (query_offset ^ subject_offset) & 1 == 0 {
        entry.set_flag(1);
        entry.set_last_hit(subject_offset.wrapping_add(diag_offset).wrapping_add(5));
    } else {
        entry.set_last_hit(subject_offset.wrapping_add(diag_offset));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-516 (the hits of a
// subject in order, each on its entry; the test of each, :518-552, is cited at `x_test_hit`)
/// LOSAT_X_BXCROSS's cell test with the kernel `Kn` (LOSAT_X_BXCROSS's AVX2 kernel, or
/// LOSAT_X_BXCROSS128's) over `cells` as for one subject's scan, its calls deferred with
/// `DEFER` (LOSAT_X_BXCROSSDEFER): every cell staged with the batches run at `CHUNK` staged
/// hits, and the stage drained after the last; checks the cells and hits counted and returns
/// the candidates and the windows with a repeated entry counted. The caller checks that the CPU
/// has the kernel's instructions ([`cross_modes`]).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn run_cells_cross<Kn: XCrossKernel, const DEFER: bool, const CHUNK: usize, F>(
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    cells: &[(i32, Vec<i32>)],
    candidate: &mut F,
) -> (u32, u64)
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    assert!(
        cells
            .iter()
            .all(|(_, hits)| hits.iter().all(|q| (0..=table.mask).contains(q))),
        "the stage recovers query offsets up to the mask only"
    );
    let mut stage = Vec::new();
    let cross = XCross::<CHUNK>::new::<DEFER>(&mut stage);
    let mut counts = XCrossCounts::default();
    let mut staged = 0usize;
    for (s, hits) in cells {
        // SAFETY: the CPU has the kernel's instructions (the caller's check), the query
        // offsets are at most the mask (asserted above), and `staged` is the count the last
        // call returned (0 at first).
        staged = unsafe {
            x_cross_cell::<Kn, i32, DEFER, CHUNK, F>(
                cross,
                staged,
                hits,
                *s,
                table,
                test,
                candidate,
                &mut counts,
            )
        };
    }
    // SAFETY: as above; `staged` is below the chunk, and the stage made for `DEFER`.
    let left = unsafe {
        Kn::drain::<DEFER, CHUNK, F>(cross, staged, true, table, test, candidate, &mut counts)
    };
    assert_eq!(left, 0, "the end of the subject drains the stage");
    assert_eq!(counts.diag.cells, cells.len() as u64, "cells staged");
    let hits: usize = cells.iter().map(|(_, hits)| hits.len()).sum();
    assert_eq!(counts.diag.hits, hits as u64, "hits staged");
    (counts.diag.candidates as u32, counts.splits)
}

// No NCBI counterpart: a test helper that picks an instance of LOSAT's stage.
/// [`run_cells_cross`] with the kernel `Kn`, its calls deferred with `DEFER`
/// (LOSAT_X_BXCROSSDEFER; the AVX2 kernel only), and the instance of `chunk`; returns the
/// candidates and the windows with a repeated entry counted.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn run_cells_cross_in<Kn: XCrossKernel, const DEFER: bool, F>(
    chunk: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    cells: &[(i32, Vec<i32>)],
    candidate: &mut F,
) -> (u32, u64)
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    match chunk {
        X_CROSS_TEST_CHUNK => {
            run_cells_cross::<Kn, DEFER, X_CROSS_TEST_CHUNK, F>(table, test, cells, candidate)
        }
        X_CROSS_CHUNK => {
            run_cells_cross::<Kn, DEFER, X_CROSS_CHUNK, F>(table, test, cells, candidate)
        }
        _ => panic!("no instance of chunk {chunk}"),
    }
}

// No NCBI counterpart: a test helper that picks LOSAT's kernel ([`run_cells_cross`] holds the
// NCBI reference).
/// [`run_cells_cross`] for a cross mode of [`cross_modes`] (LOSAT_X_BXCROSS's AVX2 kernel or
/// LOSAT_X_BXCROSS128's, its calls deferred or not) with the instance of `chunk`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn run_cells_cross_with<F>(
    mode: XCellTest,
    chunk: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    cells: &[(i32, Vec<i32>)],
    candidate: &mut F,
) -> (u32, u64)
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    match mode {
        #[cfg(target_arch = "x86_64")]
        XCellTest::Cross { defer } => {
            assert!(x_cross_available());
            if defer {
                run_cells_cross_in::<XCrossAvx2, true, F>(chunk, table, test, cells, candidate)
            } else {
                run_cells_cross_in::<XCrossAvx2, false, F>(chunk, table, test, cells, candidate)
            }
        }
        // The 128-bit kernel never defers its calls (its `drain` refuses `DEFER`).
        XCellTest::Cross128 => {
            assert!(x_cross128_available());
            run_cells_cross_in::<XCrossW128, false, F>(chunk, table, test, cells, candidate)
        }
        XCellTest::Diag(mode) => panic!("{mode:?} is not a cross mode"),
    }
}

// No NCBI counterpart: LOSAT's cell tests this machine runs.
/// The cross modes this machine runs: LOSAT_X_BXCROSS's AVX2 kernel where the CPU has AVX2
/// and POPCNT, LOSAT_X_BXCROSS128's where it has its instructions (always on wasm32 simd128),
/// each with its calls deferred (LOSAT_X_BXCROSSDEFER) and not.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn cross_modes() -> Vec<XCellTest> {
    cell_modes()
        .into_iter()
        .filter(|mode| !matches!(mode, XCellTest::Diag(_)))
        .collect()
}

/// LOSAT_X_BXCROSS's windows with a repeated entry for one subject's `cells` with the instance
/// of `chunk` ([`run_cells_cross`] on a copy of `entries`), the same for every cross mode of
/// [`cross_modes`] (both kernels, the calls deferred, LOSAT_X_BXCROSSDEFER, and not). They
/// depend on the entry coordinates only, so the test can work them out from the hits
/// ([`CrossCases`]).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn cross_splits(
    entries: &[DiagStruct],
    mask: i32,
    test: XDiagTest,
    cells: &[(i32, Vec<i32>)],
    chunk: usize,
) -> u64 {
    let splits: Vec<u64> = cross_modes()
        .into_iter()
        .map(|mode| {
            let mut entries = entries.to_vec();
            let mut table = XDiagTable::new(&mut entries, mask);
            let mut candidate = |q: i32, s: i32, _: i32, _: i32, entry: &mut DiagStruct| {
                extend(q, s, entry, test.diag_offset);
            };
            run_cells_cross_with(mode, chunk, &mut table, test, cells, &mut candidate).1
        })
        .collect();
    assert!(
        splits.windows(2).all(|pair| pair[0] == pair[1]),
        "chunk {chunk}: windows of every cross mode {splits:?}"
    );
    splits[0]
}

/// Tests the cells `(subject_offset, hits)` in order with `mode`; returns the closure's calls.
/// A batch mode's candidate count must be its number of calls. LOSAT_X_BXCROSS's mode runs
/// [`run_cells_cross`] with the instance of `chunk` (the other modes ignore it).
fn run_cells(
    mode: XCellTest,
    chunk: usize,
    entries: &mut [DiagStruct],
    mask: i32,
    test: XDiagTest,
    cells: &[(i32, Vec<i32>)],
) -> Vec<Call> {
    let mut calls = Vec::new();
    let mut counted = 0u32;
    #[cfg(target_arch = "x86_64")]
    let mut buffer = Vec::new();
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    let _ = chunk;
    let mut table = XDiagTable::new(entries, mask);
    let mut candidate = |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
        calls.push((q, s, last_hit, diff, entry.raw_bits()));
        extend(q, s, entry, test.diag_offset);
    };
    let diag_mode = match mode {
        XCellTest::Diag(mode) => Some(mode),
        #[cfg(target_arch = "x86_64")]
        XCellTest::Cross { .. } => None,
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        XCellTest::Cross128 => None,
    };
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    if diag_mode.is_none() {
        counted = run_cells_cross_with(mode, chunk, &mut table, test, cells, &mut candidate).0;
    }
    for (s, hits) in cells {
        // LOSAT_X_BXCROSS's cell test ran above.
        let Some(mode) = diag_mode else { break };
        match mode {
            XDiagMode::PerHit => x_cell_per_hit(hits, *s, &mut table, test, &mut candidate),
            XDiagMode::Scalar => {
                counted += x_cell_scalar(hits, *s, &mut table, test, &mut candidate);
            }
            #[cfg(target_arch = "x86_64")]
            XDiagMode::Avx2 => {
                assert!(x_avx2_available());
                // SAFETY: the CPU has AVX2 (asserted above).
                counted += unsafe { x_cell_avx2(hits, *s, &mut table, test, &mut candidate) };
            }
            #[cfg(target_arch = "x86_64")]
            XDiagMode::Avx2Defer => {
                assert!(x_avx2_available());
                // SAFETY: the CPU has AVX2 (asserted above).
                counted += unsafe {
                    x_cell_avx2_defer(hits, *s, &mut table, test, &mut candidate, &mut buffer)
                };
            }
            #[cfg(target_arch = "x86_64")]
            XDiagMode::Avx2DeferMasked => {
                assert!(x_avx2_available());
                // SAFETY: the CPU has AVX2 (asserted above).
                counted += unsafe {
                    x_cell_avx2_defer_masked(
                        hits,
                        *s,
                        &mut table,
                        test,
                        &mut candidate,
                        &mut buffer,
                    )
                };
            }
            #[cfg(any(
                target_arch = "x86_64",
                all(target_arch = "wasm32", target_feature = "simd128")
            ))]
            XDiagMode::W128 => {
                assert!(x_cross128_available());
                // SAFETY: on x86-64 the CPU has SSSE3, SSE4.1 and POPCNT (asserted above).
                // Without LOSAT_X_BPHITBF's test: these tables hold any entries, and the test
                // is exact only on those a scan leaves (module header); the test with it is
                // `the_branch_free_per_hit_loop_equals_the_per_hit_loop`'s.
                counted +=
                    unsafe { x_cell_w128_test(hits, *s, &mut table, test, None, &mut candidate) };
            }
        }
    }
    if mode != PER_HIT {
        assert_eq!(counted as usize, calls.len(), "{mode:?}: candidate count");
    }
    calls
}

/// Every mode against the per-hit test from the same table: the same calls in the same order
/// with the same arguments, and the same table bits. Returns the number of calls.
fn same_in_every_mode(
    entries: &[DiagStruct],
    mask: i32,
    test: XDiagTest,
    cells: &[(i32, Vec<i32>)],
) -> usize {
    let mut reference = entries.to_vec();
    let reference_calls = run_cells(PER_HIT, 0, &mut reference, mask, test, cells);
    for (mode, chunk) in cell_runs() {
        let mut got = entries.to_vec();
        let calls = run_cells(mode, chunk, &mut got, mask, test, cells);
        assert_eq!(
            calls, reference_calls,
            "{mode:?} chunk {chunk} {test:?}: candidate calls"
        );
        let at = got
            .iter()
            .zip(&reference)
            .position(|(a, b)| a.raw_bits() != b.raw_bits());
        assert!(
            at.is_none(),
            "{mode:?} chunk {chunk} {test:?}: entry {at:?}: {:#x} reference {:#x}",
            got[at.unwrap_or(0)].raw_bits(),
            reference[at.unwrap_or(0)].raw_bits()
        );
    }
    reference_calls.len()
}

/// `n` distinct query offsets from `base .. base + size` taken modulo `size` (a power of two),
/// so on distinct entries of a table of `size` entries and at most its mask, as the lookup's
/// offsets are (LOSAT_X_BXCROSS recovers them from the entry coordinate), in random order.
fn distinct_offsets(seed: &mut u64, size: usize, n: usize, base: i32) -> Vec<i32> {
    assert!(size.is_power_of_two());
    let mut all: Vec<i32> = (0..size as i32).collect();
    for i in 0..n {
        let j = i + next(seed) as usize % (size - i);
        all.swap(i, j);
    }
    all.truncate(n);
    all.into_iter()
        .map(|q| (q + base) & (size as i32 - 1))
        .collect()
}

/// Diagonal offsets with sums past 2^30 (the 31-bit field truncates) and past 2^31 - 1 (the
/// i32 sum wraps) for subject offsets up to a few thousand.
const DIAG_OFFSETS: [i32; 6] = [
    40,
    1_000,
    (1 << 30) - 600,
    (1 << 30) + 77,
    i32::MAX / 4,
    i32::MAX - 700,
];

// Every last-batch size alone (1-7) and after full batches (9-15, 17).
const HIT_COUNTS: [usize; 22] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 24, 64, 333, 1024,
];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[test]
fn batched_cells_equal_the_per_hit_test_on_random_tables() {
    let mut seed = 0x5eed_u64;
    let size = 1024usize;
    let mask = (size - 1) as i32;
    let mut calls = 0;
    for (window, word_size) in [(40, 3), (16, 2), (3, 5), (1, 0)] {
        for diag_offset in DIAG_OFFSETS {
            for (k, n) in HIT_COUNTS.into_iter().enumerate() {
                let test = XDiagTest {
                    diag_offset,
                    window,
                    word_size,
                    prefetch: [0, 3, 16][k % 3],
                };
                let s0 = (next(&mut seed) % 4000) as i32;
                // Half the entries any bit pattern, half near the cells' subject offsets
                // (flag set or not), so every outcome of the test occurs.
                let entries: Vec<DiagStruct> = (0..size)
                    .map(|_| {
                        let r = next(&mut seed);
                        if r & 1 == 0 {
                            entry_of(next(&mut seed))
                        } else {
                            let back = (next(&mut seed) % (window as u32 + 20)) as i32 - 8;
                            let last = s0.wrapping_add(diag_offset).wrapping_sub(back);
                            entry_of((last as u32 & X_LAST_HIT_MASK) | (r & 2) << 30)
                        }
                    })
                    .collect();
                // Later cells read what earlier cells (and their candidates) wrote.
                let cells: Vec<(i32, Vec<i32>)> = [0, 1, 1, 3, 10, 41]
                    .iter()
                    .map(|&step| {
                        let base = (next(&mut seed) % 3) as i32 * 7;
                        (s0 + step, distinct_offsets(&mut seed, size, n, base))
                    })
                    .collect();
                calls += same_in_every_mode(&entries, mask, test, &cells);
            }
        }
    }
    assert!(calls > 1000, "the tables must give candidates: {calls}");
}

/// The entry that gives the hit `(query_offset, s)` the outcome `case` of the test.
fn adversarial_entry(case: usize, s: i32, test: XDiagTest) -> u32 {
    let sum = s.wrapping_add(test.diag_offset);
    // last_hit for a hit without the flag at distance `diff`: s - diff + diag_offset
    let at_distance = |diff: i32| s.wrapping_sub(diff).wrapping_add(test.diag_offset);
    let flagged = |last: i32| 0x8000_0000 | (last as u32 & X_LAST_HIT_MASK);
    let clear = |last: i32| last as u32 & X_LAST_HIT_MASK;
    match case % 11 {
        // flag set: s + diag_offset above, equal to, below last_hit
        0 => flagged(sum.wrapping_add(1)),
        1 => flagged(sum),
        2 => flagged(sum.wrapping_sub(1)),
        // flag clear: diff at the window and word size edges, 0, negative, far
        3 => clear(at_distance(test.window)),
        4 => clear(at_distance(test.window - 1)),
        5 => clear(at_distance(test.word_size)),
        6 => clear(at_distance(test.word_size - 1)),
        7 => clear(at_distance(0)),
        8 => clear(at_distance(-3)),
        9 => clear(at_distance(1 << 29)),
        _ => flagged(at_distance(test.window)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[test]
fn batched_cells_equal_the_per_hit_test_at_every_edge_of_the_test() {
    let mut seed = 0xed9e_u64;
    let size = 2048usize;
    let mask = (size - 1) as i32;
    let mut calls = 0;
    for (window, word_size) in [(40, 3), (4, 3), (2, 2)] {
        for diag_offset in DIAG_OFFSETS {
            for s in [0, 5, 599, 600, 601, 699, 700, 701, 3000] {
                for (k, n) in HIT_COUNTS.into_iter().enumerate() {
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size,
                        prefetch: [16, 0, 1][k % 3],
                    };
                    let hits = distinct_offsets(&mut seed, size, n, 0);
                    let mut entries: Vec<DiagStruct> =
                        (0..size).map(|_| entry_of(next(&mut seed))).collect();
                    let rotate = next(&mut seed) as usize;
                    for (j, &q) in hits.iter().enumerate() {
                        let coord = ((q - s) & mask) as usize;
                        entries[coord] = entry_of(adversarial_entry(j + rotate, s, test));
                    }
                    // The same hits again one and window offsets later: the entries the first
                    // cell and its candidates wrote are tested again.
                    let cells = vec![(s, hits.clone()), (s + 1, hits.clone()), (s + window, hits)];
                    calls += same_in_every_mode(&entries, mask, test, &cells);
                }
            }
        }
    }
    assert!(calls > 1000, "the edges must give candidates: {calls}");
}

// ---- EXPERIMENT (LOSAT_X_BPMASKED): the masked batch against the deferred calls ----

/// blastp's default cell test (the masked batch) against the one `LOSAT_X_BPMASKED=0` runs
/// (LOSAT_X_BPDEFER's deferred calls), from the same table, on cells of every size, with
/// entries at every edge of the test and random bits, and later cells on the entries earlier
/// cells and their calls wrote: the same calls in the same order with the same arguments, and
/// the same bits in every entry of the table.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[test]
fn the_masked_batch_equals_the_deferred_calls() {
    if !x_avx2_available() {
        return;
    }
    let deferred = XCellTest::Diag(XDiagMode::Avx2Defer);
    let masked = XCellTest::Diag(XDiagMode::Avx2DeferMasked);
    let mut seed = 0xa5ced_u64;
    let size = 2048usize;
    let mask = (size - 1) as i32;
    let mut calls = 0;
    for (window, word_size) in [(40, 3), (4, 3), (2, 2)] {
        for diag_offset in DIAG_OFFSETS {
            for n in HIT_COUNTS {
                let test = XDiagTest {
                    diag_offset,
                    window,
                    word_size,
                    prefetch: 0,
                };
                let s0 = (next(&mut seed) % 4000) as i32;
                let entries: Vec<DiagStruct> = (0..size)
                    .map(|j| {
                        let r = next(&mut seed);
                        if r & 3 == 0 {
                            entry_of(next(&mut seed))
                        } else {
                            entry_of(adversarial_entry(j + r as usize, s0, test))
                        }
                    })
                    .collect();
                let cells: Vec<(i32, Vec<i32>)> = [0, 1, 2, window]
                    .iter()
                    .map(|&step| (s0 + step, distinct_offsets(&mut seed, size, n, 0)))
                    .collect();
                let mut want_table = entries.clone();
                let want = run_cells(deferred, 0, &mut want_table, mask, test, &cells);
                let mut got_table = entries;
                let got = run_cells(masked, 0, &mut got_table, mask, test, &cells);
                assert_eq!(got, want, "{n} hits {test:?}: candidate calls");
                let at = got_table
                    .iter()
                    .zip(&want_table)
                    .position(|(a, b)| a.raw_bits() != b.raw_bits());
                assert!(at.is_none(), "{n} hits {test:?}: entry {at:?}");
                calls += want.len();
            }
        }
    }
    assert!(calls > 1000, "the tables must give candidates: {calls}");
}

// ---- EXPERIMENT (LOSAT_X_BPHITBF): the branch-free test against the per-hit test ----

/// One hit by the reference (`x_test_hit`) and by `bf` then `x_test_hit` where `bf` leaves the
/// hit: the entry's bits and the result of each, and whether `bf` took the hit.
#[allow(clippy::type_complexity)]
fn hit_outcomes(
    raw: u32,
    s: i32,
    test: XDiagTest,
    bf: XHitBf,
) -> ((u32, Option<(i32, i32)>), (u32, Option<(i32, i32)>), bool) {
    let (offset, window, word) = (test.diag_offset, test.window, test.word_size);
    let mut reference = entry_of(raw);
    let expected = x_test_hit(&mut reference, s, offset, window, word);
    let mut got = entry_of(raw);
    let taken = bf.test(&mut got, bf.here(s));
    let result = if taken {
        None
    } else {
        x_test_hit(&mut got, s, offset, window, word)
    };
    (
        (reference.raw_bits(), expected),
        (got.raw_bits(), result),
        taken,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[test]
fn the_branch_free_test_gives_the_reference_outcomes_in_its_range() {
    // too far, overlap, candidate, flag 1 kept, flag 1 reset
    let mut reached = [0usize; 5];
    let pairs = [
        (2, 1),
        (3, 1),
        (3, 2),
        (40, 3),
        (41, 40),
        (1 << 20, 6),
        ((1 << 28) - 1, 1),
        ((1 << 28) - 1, (1 << 28) - 2),
    ];
    // A one-word window, a word size of 0 or of the window or more: no test.
    for (window, word_size) in [(3, 3), (40, 0), (2, 3), (1, 0)] {
        assert!(XHitBf::new(40, 300, window, word_size).is_none());
    }
    for (window, word_size) in pairs {
        // The longest subject the rule allows (at the offset 0).
        let longest = ((1i64 << 30) - 2 * i64::from(window) - 1) as usize;
        for length in [1usize, 2, 300, 1 << 20, longest] {
            // The most diagonal offset the rule allows for this subject, and no more.
            let most = ((1i64 << 30) - 2 * i64::from(window) - length as i64 - 1) as i32;
            assert!(XHitBf::new(most + 1, length, window, word_size).is_none());
            assert!(XHitBf::new(-1, length, window, word_size).is_none());
            let offsets = [0, 1, window, most / 2, most];
            for diag_offset in offsets.into_iter().filter(|&offset| offset <= most) {
                let bf = XHitBf::new(diag_offset, length, window, word_size).expect("in range");
                let test = XDiagTest {
                    diag_offset,
                    window,
                    word_size,
                    prefetch: 0,
                };
                // The flag-0 values a table holds (module header): [-window, diag_offset +
                // length).
                let (lowest, highest) = (-window, diag_offset + length as i32 - 1);
                for s in [0, length as i32 - 1] {
                    let here = s + diag_offset;
                    let flag0 = [
                        here - (word_size - 1),
                        here - word_size,
                        here - (window - 1),
                        here - window,
                        here - window - 1,
                        here,
                        here + 1,
                        0,
                        lowest,
                        highest,
                    ];
                    let flag1 = [here + 1, here, here - 1, lowest, highest];
                    let entries = flag0
                        .into_iter()
                        .filter(|last| (lowest..=highest).contains(last))
                        .map(|last| last as u32 & X_LAST_HIT_MASK)
                        .chain(flag1.map(|last| 0x8000_0000 | (last as u32 & X_LAST_HIT_MASK)));
                    for raw in entries {
                        let (expected, got, taken) = hit_outcomes(raw, s, test, bf);
                        assert_eq!(
                            got, expected,
                            "window {window} word {word_size} length {length} offset {diag_offset} s {s} raw {raw:#x}"
                        );
                        let outcome = match (raw >> 31, expected) {
                            (0, (_, Some(_))) => 2,
                            (0, (after, None)) if after != raw => 0,
                            (0, _) => 1,
                            (_, (after, _)) if after == raw => 3,
                            _ => 4,
                        };
                        assert_eq!(taken, outcome < 2, "raw {raw:#x}: taken");
                        reached[outcome] += 1;
                    }
                }
            }
        }
    }
    assert!(reached.iter().all(|&n| n > 0), "every outcome: {reached:?}");
    // Just past the range the doubled difference wraps: a cleared entry at the offset 2^30 is
    // too far for the reference and an overlap for the test.
    let (window, word_size, diag_offset) = (40, 3, 1 << 30);
    let bf = XHitBf {
        diag_offset,
        word2: 2 * word_size,
        span2: 2 * (window - word_size) as u32,
        window2: 2 * window,
    };
    let test = XDiagTest {
        diag_offset,
        window,
        word_size,
        prefetch: 0,
    };
    let (expected, got, _) = hit_outcomes(-window as u32 & X_LAST_HIT_MASK, 0, test, bf);
    assert_ne!(got, expected, "outside the range");
}

/// [`run_cells`] with the per-hit loop of the branch-free test.
fn run_cells_bf(
    entries: &mut [DiagStruct],
    mask: i32,
    test: XDiagTest,
    bf: XHitBf,
    cells: &[(i32, Vec<i32>)],
) -> Vec<Call> {
    let mut calls = Vec::new();
    let mut table = XDiagTable::new(entries, mask);
    let mut candidate = |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
        calls.push((q, s, last_hit, diff, entry.raw_bits()));
        extend(q, s, entry, test.diag_offset);
    };
    for (s, hits) in cells {
        x_cell_per_hit_bf(hits, *s, &mut table, test, bf, &mut candidate);
    }
    calls
}

/// LOSAT_X_BPCELL128's cell test with LOSAT_X_BPHITBF's test for the last hits
/// ([`x_cell_w128`] with `bf`) over `cells`; checks its count of candidates.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn run_cells_w128_bf(
    entries: &mut [DiagStruct],
    mask: i32,
    test: XDiagTest,
    bf: XHitBf,
    cells: &[(i32, Vec<i32>)],
) -> Vec<Call> {
    assert!(x_cross128_available());
    let mut calls = Vec::new();
    let mut counted = 0u32;
    let mut table = XDiagTable::new(entries, mask);
    let mut candidate = |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
        calls.push((q, s, last_hit, diff, entry.raw_bits()));
        extend(q, s, entry, test.diag_offset);
    };
    for (s, hits) in cells {
        // SAFETY: on x86-64 the CPU has SSSE3, SSE4.1 and POPCNT (asserted above); the
        // offsets are below the table's size, and the caller's `bf` fits the subject.
        counted +=
            unsafe { x_cell_w128_test(hits, *s, &mut table, test, Some(bf), &mut candidate) };
    }
    assert_eq!(counted as usize, calls.len(), "{test:?}: candidate count");
    calls
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[test]
fn the_branch_free_per_hit_loop_equals_the_per_hit_loop() {
    let mut seed = 0xb7_u64;
    let size = 1024usize;
    let mask = (size - 1) as i32;
    let length = 5_000usize;
    let mut calls = 0;
    for (window, word_size) in [(40, 3), (16, 2), (4, 3), (2, 1)] {
        // From the first offset to the most the rule allows for the subject.
        let most = ((1i64 << 30) - 2 * i64::from(window) - length as i64 - 1) as i32;
        for diag_offset in [window, 1_000, i32::MAX / 4, most] {
            let bf = XHitBf::new(diag_offset, length, window, word_size).expect("in range");
            let (lowest, highest) = (-window, diag_offset + length as i32 - 1);
            for (k, n) in HIT_COUNTS.into_iter().enumerate() {
                // The read-ahead distances of the other tests, and one equal to the cell's hits
                // and one past them (no hit reads ahead).
                let test = XDiagTest {
                    diag_offset,
                    window,
                    word_size,
                    prefetch: [0, 3, 16, n, n + 1][k % 5],
                };
                let s0 = (next(&mut seed) % 4000) as i32;
                // The entries a table holds (module header): flag 0 with last_hit in [-window,
                // diag_offset + length), most of them near the cells' subject offsets, or
                // flag 1.
                let entries: Vec<DiagStruct> = (0..size)
                    .map(|_| {
                        let r = next(&mut seed);
                        let last = if r.is_multiple_of(4) {
                            lowest + (next(&mut seed) % (highest - lowest + 1) as u32) as i32
                        } else {
                            let back = (next(&mut seed) % (window as u32 + 20)) as i32 - 8;
                            (s0 + diag_offset - back).clamp(lowest, highest)
                        };
                        let flag = if r % 4 == 1 { 0x8000_0000 } else { 0 };
                        entry_of(flag | (last as u32 & X_LAST_HIT_MASK))
                    })
                    .collect();
                // Later cells read what earlier cells (and their candidates) wrote.
                let cells: Vec<(i32, Vec<i32>)> = [0, 1, 1, 3, 10, 41]
                    .iter()
                    .map(|&step| {
                        let base = (next(&mut seed) % 3) as i32 * 7;
                        (s0 + step, distinct_offsets(&mut seed, size, n, base))
                    })
                    .collect();
                let mut reference = entries.clone();
                let expected = run_cells(PER_HIT, 0, &mut reference, mask, test, &cells);
                let mut got = entries.clone();
                let got_calls = run_cells_bf(&mut got, mask, test, bf, &cells);
                assert_eq!(got_calls, expected, "{test:?}: candidate calls");
                let at = got
                    .iter()
                    .zip(&reference)
                    .position(|(a, b)| a.raw_bits() != b.raw_bits());
                assert!(at.is_none(), "{test:?}: entry {at:?}");
                // LOSAT_X_BPCELL128: the batches of eight, then the same test for the last
                // hits of each cell.
                #[cfg(any(
                    target_arch = "x86_64",
                    all(target_arch = "wasm32", target_feature = "simd128")
                ))]
                if x_cross128_available() {
                    let mut got = entries.clone();
                    let got_calls = run_cells_w128_bf(&mut got, mask, test, bf, &cells);
                    assert_eq!(got_calls, expected, "{test:?}: 128-bit candidate calls");
                    let at = got
                        .iter()
                        .zip(&reference)
                        .position(|(a, b)| a.raw_bits() != b.raw_bits());
                    assert!(at.is_none(), "{test:?}: 128-bit entry {at:?}");
                }
                calls += expected.len();
            }
        }
    }
    assert!(calls > 1000, "the tables must give candidates: {calls}");
}

/// One subject's scan with the cell test `mode` and the scan of `branch_free`
/// ([`x_scan_two_hit_in_with`]; BXSCANBF, graduated for blastx), LOSAT_X_BXCROSS's cell test as blastx's
/// entry runs it ([`x_scan_cells_cross`]) with the instance of `chunk`; returns the counts of
/// the cells tested in batches.
#[allow(clippy::too_many_arguments)]
fn scan_subject<K: XCell, F>(
    (mode, chunk, branch_free): (XCellTest, usize, bool),
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> XDiagCounts
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    #[cfg_attr(
        not(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        )),
        expect(
            clippy::infallible_destructuring_match,
            reason = "the `Cross` and `Cross128` variants exist on x86_64 and wasm32 simd128 only, so elsewhere `Diag` is the only pattern"
        )
    )]
    let mode = match mode {
        XCellTest::Diag(mode) => mode,
        #[cfg(target_arch = "x86_64")]
        XCellTest::Cross { defer } => {
            return match chunk {
                X_CROSS_TEST_CHUNK => x_scan_cells_cross::<X_CROSS_TEST_CHUNK, K, F>(
                    branch_free,
                    defer,
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                ),
                X_CROSS_CHUNK => x_scan_cells_cross::<X_CROSS_CHUNK, K, F>(
                    branch_free,
                    defer,
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                ),
                _ => panic!("no instance of chunk {chunk}"),
            }
            .diag;
        }
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        XCellTest::Cross128 => {
            return match chunk {
                X_CROSS_TEST_CHUNK => x_scan_cells_cross128::<X_CROSS_TEST_CHUNK, K, F>(
                    branch_free,
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                ),
                X_CROSS_CHUNK => x_scan_cells_cross128::<X_CROSS_CHUNK, K, F>(
                    branch_free,
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                ),
                _ => panic!("no instance of chunk {chunk}"),
            }
            .diag;
        }
    };
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    let _ = chunk;
    if branch_free {
        x_scan_two_hit_in_with::<true, _, F>(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            candidate,
        )
    } else {
        x_scan_two_hit_in_with::<false, _, F>(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            candidate,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-512 and c++/src/algo/blast/core/blast_aascan.c:60-99
// (cited at `x_scan_two_hit` and `x_scan_cells`)
#[test]
fn batched_scans_of_real_lookup_tables_equal_the_per_hit_scan() {
    use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
    use crate::algorithm::tblastx::translation::QueryFrame;
    use crate::config::ScoringMatrix;
    use crate::stats::lookup_protein_params_ungapped;

    // NCBISTDAA codes of the 20 standard amino acids
    const STANDARD: [u8; 20] = [
        1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
    ];
    let mut seed = 0x1007_u64;
    let mut protein = |n: usize| -> Vec<u8> {
        (0..n)
            .map(|_| STANDARD[next(&mut seed) as usize % STANDARD.len()])
            .collect()
    };
    let mut queries: Vec<Vec<u8>> = (0..6).map(|i| protein(80 + 90 * i)).collect();
    // Low complexity: cells far above AA_HITS_PER_CELL (the overflow array).
    let mut low = vec![11u8; 150];
    low.extend(protein(30));
    low.extend(std::iter::repeat_n([17u8, 18], 60).flatten());
    queries.push(low);
    let frames: Vec<Vec<QueryFrame>> = queries
        .iter()
        .map(|seq| {
            let mut aa_seq = vec![0u8];
            aa_seq.extend_from_slice(seq);
            aa_seq.push(0);
            vec![QueryFrame {
                frame: 1,
                aa_seq,
                aa_seq_nomask: None,
                aa_len: seq.len(),
                orig_len: seq.len(),
                seg_masks: Vec::new(),
            }]
        })
        .collect();
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (lookup, _) = build_ncbi_lookup(&frames, 11, &params, false);
    assert!(
        lookup
            .backbone
            .iter()
            .any(|cell| cell.num_used as usize > 3 * X_DIAG_LANES),
        "the query set must give long cells"
    );
    // Subjects: unrelated proteins, copies of the queries with every 5th letter changed (two
    // hits on a diagonal, extensions), the low-complexity query, and short ones.
    let mut subjects: Vec<Vec<u8>> = queries
        .iter()
        .map(|q| {
            q.iter()
                .enumerate()
                .map(|(i, &a)| if i % 5 == 4 { STANDARD[i % 20] } else { a })
                .collect()
        })
        .collect();
    subjects.extend((0..4).map(|i| protein(50 + 200 * i)));
    subjects.push(queries[6].clone());
    subjects.push(protein(1));
    subjects.push(protein(3));
    let window = 40;
    let mut size = 1usize;
    while size < lookup.query_length as usize + window as usize {
        size <<= 1;
    }
    let mask = (size - 1) as i32;
    let mut calls = 0;
    for start in [window, (1 << 30) - 2_000, i32::MAX - 3_000] {
        for prefetch in [16, 0, 5] {
            let mut results = Vec::new();
            // Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx): every cell test also with the branch-free scan
            // (LOSAT_X_BXCROSS's with each chunk).
            let runs = cell_runs()
                .into_iter()
                .flat_map(|(mode, chunk)| [(mode, chunk, false), (mode, chunk, true)]);
            for (mode, chunk, branch_free) in runs {
                let mut entries = vec![DiagStruct::clear(window); size];
                let mut log: Vec<(usize, Call)> = Vec::new();
                let mut diag_offset: i32 = start;
                let mut counts = (0u64, 0u64, 0u64);
                for (k, subject) in subjects.iter().enumerate() {
                    let mut buffer = subject.clone();
                    buffer.push(0);
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size: 3,
                        prefetch,
                    };
                    let mut candidate =
                        |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
                            log.push((k, (q, s, last_hit, diff, entry.raw_bits())));
                            extend(q, s, entry, diag_offset);
                        };
                    let table = &mut XDiagTable::new(&mut entries, mask);
                    let c = scan_subject(
                        (mode, chunk, branch_free),
                        &lookup,
                        &buffer,
                        subject.len(),
                        table,
                        test,
                        &mut candidate,
                    );
                    counts = (
                        counts.0 + c.cells,
                        counts.1 + c.hits,
                        counts.2 + c.candidates,
                    );
                    diag_offset = diag_offset.wrapping_add(subject.len() as i32 + window);
                }
                if mode == PER_HIT {
                    assert_eq!(counts, (0, 0, 0), "the per-hit loop counts nothing");
                } else {
                    assert!(counts.0 > 0 && counts.1 >= counts.0, "{mode:?}: {counts:?}");
                    assert_eq!(counts.2 as usize, log.len(), "{mode:?}: candidates");
                }
                let bits: Vec<u32> = entries.iter().map(|e| e.raw_bits()).collect();
                results.push(((mode, chunk, branch_free), log, bits, counts));
            }
            let (_, reference_log, reference_bits, _) = &results[0];
            // Every batch mode counts the same cells, hits and candidates (LOSAT_X_BXCROSS's
            // as it stages the cells and drains the stage, the others per cell).
            let batched = results
                .iter()
                .find(|(run, ..)| run.0 != PER_HIT)
                .map(|(.., counts)| *counts);
            for (mode, log, bits, counts) in &results[1..] {
                assert!(
                    log == reference_log,
                    "{mode:?} start {start} prefetch {prefetch}: candidate calls differ"
                );
                assert!(
                    bits == reference_bits,
                    "{mode:?} start {start} prefetch {prefetch}: tables differ"
                );
                if mode.0 != PER_HIT {
                    assert_eq!(Some(*counts), batched, "{mode:?}: counts");
                }
            }
            calls += reference_log.len();
        }
    }
    assert!(calls > 100, "the subjects must give candidates: {calls}");
}

/// The query set of the tests of the table of 16-bit hits, its tables of 32 and 16 bits (cells
/// of every length: eight-lane tails of 1-7 hits, cells in the overflow array and the padding
/// at its end), the subjects to scan and the diagonal table's `(size, mask, window)`.
#[allow(clippy::type_complexity)]
fn tables_of_both_widths() -> (
    BlastAaLookupTable,
    BlastAaLookupTable<crate::algorithm::tblastx::lookup::BackboneCell16>,
    Vec<Vec<u8>>,
    (usize, i32, i32),
) {
    use crate::algorithm::tblastx::lookup::XBpLookup;
    use crate::algorithm::tblastx::lookup::{x_bp_build_ncbi_lookup, x_bp_build_ncbi_lookup_as};
    use crate::algorithm::tblastx::translation::QueryFrame;

    const STANDARD: [u8; 20] = [
        1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
    ];
    let mut seed = 0x16_u64;
    let mut protein = |n: usize| -> Vec<u8> {
        (0..n)
            .map(|_| STANDARD[next(&mut seed) as usize % STANDARD.len()])
            .collect()
    };
    let mut queries: Vec<Vec<u8>> = (0..5).map(|i| protein(70 + 110 * i)).collect();
    let mut low = vec![11u8; 150];
    low.extend(protein(30));
    low.extend(std::iter::repeat_n([17u8, 18], 60).flatten());
    queries.push(low);
    let frames: Vec<Vec<QueryFrame>> = queries
        .iter()
        .map(|seq| {
            let mut aa_seq = vec![0u8];
            aa_seq.extend_from_slice(seq);
            aa_seq.push(0);
            vec![QueryFrame {
                frame: 1,
                aa_seq,
                aa_seq_nomask: None,
                aa_len: seq.len(),
                orig_len: seq.len(),
                seg_masks: Vec::new(),
            }]
        })
        .collect();
    let (wide, _) = x_bp_build_ncbi_lookup(&frames, 11, false);
    let (XBpLookup::Narrow(narrow), _) = x_bp_build_ncbi_lookup_as(
        &frames,
        11,
        false,
        true,
        crate::config::ScoringMatrix::Blosum62,
        crate::algorithm::tblastx::lookup::LOOKUP_WORD_LENGTH,
    ) else {
        panic!("the query set fits 16 bits");
    };
    let lengths: Vec<usize> = narrow.backbone.iter().map(|cell| cell.count()).collect();
    for tail in 1..X_DIAG_LANES {
        assert!(
            lengths.iter().any(|&n| n > 3 && n % X_DIAG_LANES == tail),
            "the query set must give a cell with {tail} hits past a batch of eight"
        );
    }
    assert!(lengths.iter().any(|&n| n > 3 * X_DIAG_LANES));
    let mut subjects: Vec<Vec<u8>> = queries
        .iter()
        .map(|q| {
            q.iter()
                .enumerate()
                .map(|(i, &a)| if i % 5 == 4 { STANDARD[i % 20] } else { a })
                .collect()
        })
        .collect();
    subjects.extend((0..3).map(|i| protein(60 + 300 * i)));
    subjects.push(queries[5].clone());
    subjects.push(protein(1));
    subjects.push(protein(3));
    let window = 40;
    let mut size = 1usize;
    while size < wide.query_length as usize + window as usize {
        size <<= 1;
    }
    let mask = (size - 1) as i32;
    (wide, narrow, subjects, (size, mask, window))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-512 and c++/src/algo/blast/core/blast_aascan.c:60-99
// (cited at `x_scan_two_hit` and `x_scan_cells`)
// LOSAT_X_BPLK16: the scans of the table of 16-bit hits (every cell test of blastp's loops, with
// both scans of the cells, and the scans that list the cells) call `candidate` and `on_cell`
// with the same offsets in the same order, and leave the same diagonal entries, as the scans of
// the table of 32-bit hits of the same query, on a query set with cells of every length
// (eight-lane tails of 1-7 hits, cells in the overflow array and the padding at its end).
#[test]
fn scans_of_the_table_of_u16_hits_equal_the_scans_of_the_table_of_i32_hits() {
    let (wide, narrow, subjects, (size, mask, window)) = tables_of_both_widths();
    // The calls, the entries and the cells of every cell test and scan over the subjects.
    fn run<K: XCell>(
        lookup: &BlastAaLookupTable<K>,
        subjects: &[Vec<u8>],
        (size, mask, window): (usize, i32, i32),
    ) -> Vec<(Vec<(usize, Call)>, Vec<u32>)> {
        let mut runs = Vec::new();
        for mode in cell_modes() {
            #[cfg_attr(
                not(any(
                    target_arch = "x86_64",
                    all(target_arch = "wasm32", target_feature = "simd128")
                )),
                expect(
                    irrefutable_let_patterns,
                    reason = "the `Cross` and `Cross128` variants exist on x86_64 and wasm32 simd128 only, so elsewhere `Diag` is the only pattern"
                )
            )]
            let XCellTest::Diag(mode) = mode
            else {
                continue;
            };
            for branch_free in [false, true] {
                let mut entries = vec![DiagStruct::clear(window); size];
                let mut log: Vec<(usize, Call)> = Vec::new();
                let mut diag_offset: i32 = window;
                for (k, subject) in subjects.iter().enumerate() {
                    let mut buffer = subject.clone();
                    buffer.push(0);
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size: 3,
                        prefetch: 16,
                    };
                    let mut candidate =
                        |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
                            log.push((k, (q, s, last_hit, diff, entry.raw_bits())));
                            extend(q, s, entry, diag_offset);
                        };
                    let table = &mut XDiagTable::new(&mut entries, mask);
                    if branch_free {
                        x_scan_two_hit_in_with::<true, _, _>(
                            mode,
                            lookup,
                            &buffer,
                            subject.len(),
                            table,
                            test,
                            &mut candidate,
                        );
                    } else {
                        x_scan_two_hit_in_with::<false, _, _>(
                            mode,
                            lookup,
                            &buffer,
                            subject.len(),
                            table,
                            test,
                            &mut candidate,
                        );
                    }
                    diag_offset = diag_offset.wrapping_add(subject.len() as i32 + window);
                }
                runs.push((log, entries.iter().map(|e| e.raw_bits()).collect()));
            }
        }
        // The per-hit loop with the branch-free two-hit test (LOSAT_X_BPHITBF).
        let mut entries = vec![DiagStruct::clear(window); size];
        let mut log: Vec<(usize, Call)> = Vec::new();
        let mut diag_offset: i32 = window;
        for (k, subject) in subjects.iter().enumerate() {
            let mut buffer = subject.clone();
            buffer.push(0);
            let test = XDiagTest {
                diag_offset,
                window,
                word_size: 3,
                prefetch: 16,
            };
            let bf = XHitBf::new(diag_offset, subject.len(), window, 3).expect("in range");
            let mut candidate =
                |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
                    log.push((k, (q, s, last_hit, diff, entry.raw_bits())));
                    extend(q, s, entry, diag_offset);
                };
            let table = &mut XDiagTable::new(&mut entries, mask);
            x_scan_two_hit_hitbf(
                lookup,
                &buffer,
                subject.len(),
                table,
                test,
                bf,
                &mut candidate,
            );
            diag_offset = diag_offset.wrapping_add(subject.len() as i32 + window);
        }
        runs.push((log, entries.iter().map(|e| e.raw_bits()).collect()));
        runs
    }
    let want = run(&wide, &subjects, (size, mask, window));
    let got = run(&narrow, &subjects, (size, mask, window));
    assert_eq!(want.len(), got.len());
    let mut calls = 0;
    for (i, (want, got)) in want.iter().zip(&got).enumerate() {
        assert!(want.0 == got.0, "run {i}: candidate calls differ");
        assert!(want.1 == got.1, "run {i}: diagonal entries differ");
        calls += want.0.len();
    }
    assert!(calls > 100, "the subjects must give candidates: {calls}");
    // The scans that list the cells.
    let mut compared = 0usize;
    for subject in &subjects {
        let mut buffer = subject.clone();
        buffer.push(0);
        let (mut a, mut b, mut c): (Vec<(i32, Vec<i32>)>, Vec<_>, Vec<_>) =
            (Vec::new(), Vec::new(), Vec::new());
        x_scan_cells(&wide, &buffer, subject.len(), |s, hits| {
            a.push((s, hits.to_vec()))
        });
        x_scan_cells(&narrow, &buffer, subject.len(), |s, hits| {
            b.push((s, hits.iter().map(|h| h.q()).collect()))
        });
        x_scan_cells_with::<true, _, _>(&narrow, &buffer, subject.len(), |s, hits| {
            c.push((s, hits.iter().map(|h| h.q()).collect()))
        });
        assert!(
            a == b && a == c,
            "the cells of a subject of {}",
            subject.len()
        );
        compared += a.len();
    }
    assert!(compared > 100, "{compared} cells compared");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-512 and c++/src/algo/blast/core/blast_aascan.c:60-99
// (cited at `x_scan_two_hit` and `x_scan_cells`)
// BLASTX's table of 16-bit hits (BackboneCell16): every cell test of blastx's fused loop,
// LOSAT_X_BXCROSS's among them (with the instances of both chunks, its calls deferred and not,
// and both scans of the cells), calls `candidate` with the same arguments in the same order,
// counts the same cells, hits and candidates and windows with a repeated entry, and leaves the
// same diagonal entries on the table of 16 bits as on the table of 32 bits of the same query, and
// as the per-hit loop on that table. The query set has cells of 1-7 hits past a batch of eight
// and cells in the overflow array (LOSAT_X_BXCROSS's three loads of eight hits read 48 bytes
// from the start of a list, into the padding at the end of the cells and of the overflow).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[test]
fn the_cross_cell_test_on_u16_hits_equals_the_cell_tests_on_i32_hits() {
    let (wide, narrow, subjects, (size, mask, window)) = tables_of_both_widths();
    // The reference: the per-hit loop on the table of 32 bits.
    fn run<K: XCell>(
        lookup: &BlastAaLookupTable<K>,
        subjects: &[Vec<u8>],
        (size, mask, window): (usize, i32, i32),
        run: (XCellTest, usize, bool),
    ) -> (Vec<(usize, Call)>, Vec<u32>, [u64; 4]) {
        let mut entries = vec![DiagStruct::clear(window); size];
        let mut log: Vec<(usize, Call)> = Vec::new();
        let mut diag_offset: i32 = window;
        let mut counts = [0u64; 4];
        for (k, subject) in subjects.iter().enumerate() {
            let mut buffer = subject.clone();
            buffer.push(0);
            let test = XDiagTest {
                diag_offset,
                window,
                word_size: 3,
                prefetch: 16,
            };
            let mut candidate =
                |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
                    log.push((k, (q, s, last_hit, diff, entry.raw_bits())));
                    extend(q, s, entry, diag_offset);
                };
            let table = &mut XDiagTable::new(&mut entries, mask);
            let c = scan_subject(
                run,
                lookup,
                &buffer,
                subject.len(),
                table,
                test,
                &mut candidate,
            );
            for (total, n) in counts.iter_mut().zip([c.cells, c.hits, c.candidates, 0]) {
                *total += n;
            }
            diag_offset = diag_offset.wrapping_add(subject.len() as i32 + window);
        }
        (log, entries.iter().map(|e| e.raw_bits()).collect(), counts)
    }
    let reference = run(&wide, &subjects, (size, mask, window), (PER_HIT, 0, false));
    assert!(reference.0.len() > 100, "{} candidates", reference.0.len());
    #[cfg(target_arch = "x86_64")]
    let mut crossed = 0;
    let mut crossed128 = 0;
    for (mode, chunk) in cell_runs() {
        for branch_free in [false, true] {
            let what = format!("{mode:?} chunk {chunk} branch_free {branch_free}");
            let r = (mode, chunk, branch_free);
            let on_wide = run(&wide, &subjects, (size, mask, window), r);
            let on_narrow = run(&narrow, &subjects, (size, mask, window), r);
            assert!(on_narrow.0 == reference.0, "{what}: candidate calls");
            assert!(on_narrow.1 == reference.1, "{what}: diagonal entries");
            assert_eq!(on_narrow.2[..3], on_wide.2[..3], "{what}: counts");
            #[cfg(target_arch = "x86_64")]
            {
                crossed += usize::from(matches!(mode, XCellTest::Cross { .. }));
            }
            crossed128 += usize::from(mode == XCellTest::Cross128);
        }
    }
    #[cfg(target_arch = "x86_64")]
    assert!(
        crossed > 0 || !x_cross_available(),
        "the cross cell test must run where the CPU has it"
    );
    assert!(
        crossed128 > 0 || !x_cross128_available(),
        "the 128-bit cross cell test must run where the CPU has it"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:60-99 (cited at `x_scan_cells`)
// Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx): the branch-free scan calls `on_cell` with the positions and
// hit lists of the presence-test loop, on subjects across the block boundaries, with the
// letters of translated frames (X, stop; all of 1-27) and with the word sizes 2-4 (the general
// loop).
#[test]
fn the_branch_free_scan_finds_the_cells_of_the_presence_test() {
    use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
    use crate::algorithm::tblastx::translation::QueryFrame;
    use crate::config::ScoringMatrix;
    use crate::stats::lookup_protein_params_ungapped;

    let mut seed = 0x5ca7_u64;
    let mut letters = |n: usize, top: u64| -> Vec<u8> {
        (0..n)
            .map(|_| 1 + (u64::from(next(&mut seed)) % top) as u8)
            .collect()
    };
    let mut compared = 0usize;
    for (case, word_length) in [3, 3, 2, 4, 3].into_iter().enumerate() {
        let frames: Vec<Vec<QueryFrame>> = (0..1 + 3 * case)
            .map(|i| {
                let seq = letters(60 + 70 * i, 22);
                let mut aa_seq = vec![0u8];
                aa_seq.extend_from_slice(&seq);
                aa_seq.push(0);
                vec![QueryFrame {
                    frame: 1,
                    aa_seq,
                    aa_seq_nomask: None,
                    aa_len: seq.len(),
                    orig_len: seq.len(),
                    seg_masks: Vec::new(),
                }]
            })
            .collect();
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (lookup, _) = build_ncbi_lookup_for_profile(
            &frames,
            [11, 13][case % 2],
            &params,
            false,
            ScoringMatrix::Blosum62,
            word_length,
        );
        for length in [0, 1, 2, 3, 4, 255, 256, 257, 258, 511, 512, 513, 2_000] {
            let mut subject = letters(length, 27);
            subject.push(0);
            let mut want: Vec<(i32, Vec<i32>)> = Vec::new();
            let mut got: Vec<(i32, Vec<i32>)> = Vec::new();
            x_scan_cells(&lookup, &subject, length, |s, hits| {
                want.push((s, hits.to_vec()))
            });
            x_scan_cells_with::<true, _, _>(&lookup, &subject, length, |s, hits| {
                got.push((s, hits.to_vec()))
            });
            assert_eq!(got, want, "case {case} length {length}");
            compared += want.len();
        }
    }
    assert!(compared > 1_000, "{compared} cells compared");
}

#[test]
fn the_pack_table_lists_the_set_lanes_in_order() {
    for (mask, order) in X_PACK_LANES.iter().enumerate() {
        let lanes: Vec<u8> = (0..X_DIAG_LANES as u8)
            .filter(|&lane| mask >> lane & 1 != 0)
            .collect();
        assert_eq!(&order[..lanes.len()], &lanes[..], "mask {mask:#x}");
        assert!(order[lanes.len()..].iter().all(|&lane| lane == 0));
    }
}

#[test]
fn the_batch_check_sees_a_repeated_entry() {
    assert!(x_batch_coords_distinct(&[]));
    assert!(x_batch_coords_distinct(&[3, 1, 2, 0, 7, 6, 5, 4]));
    assert!(!x_batch_coords_distinct(&[3, 1, 2, 0, 7, 6, 5, 3]));
}

// ---- EXPERIMENT (LOSAT_X_BXCROSS): the batches across the cells' borders ----

/// What the cells of subjects give LOSAT_X_BXCROSS's stage, worked out from the hits and the
/// reference's calls by following the pieces and drains of `x_cross_cell` and `x_cross_drain`
/// and the parts of `x_cross_window`, so that a test can require each case to occur.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[derive(Debug, Default)]
struct CrossCases {
    /// Windows of eight with a repeated entry (tested in parts by the AVX2 kernel).
    windows: usize,
    /// Parts of such windows cut short before a repeated entry, by the hits they take (1-7).
    cuts: [usize; X_DIAG_LANES],
    /// Parts that start past their window's first lane (after a cut).
    late_parts: usize,
    /// Windows of eight distinct entries with a candidate in a lane before a lane of a later
    /// cell.
    candidate_before_later_cell: usize,
    /// Cells of more than 24 hits (more than one piece), and of more than 48.
    long_cells: usize,
    longer_cells: usize,
    /// Drains while staging (the chunk reached), and those in the middle of a cell.
    drains: usize,
    drains_in_a_cell: usize,
    /// The hits a drain while staging moves to the front of the stage, by count (0-7).
    moved: [usize; X_DIAG_LANES],
    /// The hits left for the one-by-one test at the end of a subject, by count (0-7).
    left: [usize; X_DIAG_LANES],
}

#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
impl CrossCases {
    /// Adds the cases of one subject's `cells`, whose reference calls are `calls`, with the
    /// batches run at `chunk` staged hits.
    fn add(&mut self, cells: &[(i32, Vec<i32>)], mask: i32, calls: &[Call], chunk: usize) {
        // The stream of staged hits (q, s) and their entries.
        let stream: Vec<(i32, i32)> = cells
            .iter()
            .flat_map(|(s, hits)| hits.iter().map(move |&q| (q, *s)))
            .collect();
        let coords: Vec<usize> = stream
            .iter()
            .map(|&(q, s)| (q.wrapping_sub(s) & mask) as usize)
            .collect();
        let called: std::collections::HashSet<(i32, i32)> =
            calls.iter().map(|&(q, s, ..)| (q, s)).collect();
        // x_cross_drain over the stream's hits `front..front + staged` (`front` at the stage's
        // first position): windows of eight from the front; returns the hits left.
        let drain = |cases: &mut Self, front: usize, staged: usize| {
            let full = staged - staged % X_DIAG_LANES;
            for b in (front..front + full).step_by(X_DIAG_LANES) {
                let lanes = &coords[b..b + X_DIAG_LANES];
                if (1..X_DIAG_LANES).any(|j| lanes[..j].contains(&lanes[j])) {
                    cases.windows += 1;
                    let mut start = 0;
                    while start < X_DIAG_LANES {
                        let mut seen = std::collections::HashSet::new();
                        let end = (start..X_DIAG_LANES)
                            .find(|&j| !seen.insert(lanes[j]))
                            .unwrap_or(X_DIAG_LANES);
                        if end < X_DIAG_LANES {
                            cases.cuts[end - start] += 1;
                        }
                        cases.late_parts += usize::from(start > 0);
                        start = end;
                    }
                    continue;
                }
                let batch = &stream[b..b + X_DIAG_LANES];
                if batch.iter().enumerate().any(|(i, hit)| {
                    called.contains(hit) && batch[i + 1..].iter().any(|later| later.1 != hit.1)
                }) {
                    cases.candidate_before_later_cell += 1;
                }
            }
            staged - full
        };
        let (mut front, mut staged) = (0usize, 0usize);
        for (_, hits) in cells {
            self.long_cells += usize::from(hits.len() > X_CROSS_PIECE);
            self.longer_cells += usize::from(hits.len() > 2 * X_CROSS_PIECE);
            let mut left = hits.len();
            loop {
                let piece = left.min(X_CROSS_PIECE);
                staged += piece;
                if staged >= chunk {
                    self.drains += 1;
                    self.drains_in_a_cell += usize::from(left > X_CROSS_PIECE);
                    let rest = drain(self, front, staged);
                    self.moved[rest] += 1;
                    front += staged - rest;
                    staged = rest;
                }
                if left <= X_CROSS_PIECE {
                    break;
                }
                left -= X_CROSS_PIECE;
            }
        }
        let rest = drain(self, front, staged);
        self.left[rest] += 1;
    }
}

/// The reference's calls for `cells` from `entries` (the per-hit test).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn reference_calls(
    entries: &[DiagStruct],
    mask: i32,
    test: XDiagTest,
    cells: &[(i32, Vec<i32>)],
) -> Vec<Call> {
    run_cells(PER_HIT, 0, &mut entries.to_vec(), mask, test, cells)
}

/// Entries near `s0`'s diagonal values (flag set or not) for half the table and any bits for
/// the other half, so that every outcome of the test occurs.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn entries_near(seed: &mut u64, size: usize, s0: i32, test: XDiagTest) -> Vec<DiagStruct> {
    (0..size)
        .map(|_| {
            let r = next(seed);
            if r & 1 == 0 {
                entry_of(next(seed))
            } else {
                let back = (next(seed) % (test.window as u32 + 20)) as i32 - 8;
                let last = s0.wrapping_add(test.diag_offset).wrapping_sub(back);
                entry_of((last as u32 & X_LAST_HIT_MASK) | (r & 2) << 30)
            }
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
// Subjects of 7 to 200 cells at consecutive subject offsets (some gaps), each cell of 0 to 40
// hits or of 25 to 120 (more than one piece of the copy), one of 333 in some: the batches cross
// the cells' borders, the stage is drained between cells and inside long ones (with the chunk
// of the test hook and the built-in one), every count of hits is moved to the front by a drain
// and left for the end. Later cells hit the entries earlier cells and their candidates wrote.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[test]
fn cross_batches_equal_the_per_hit_test_on_cells_of_every_size() {
    if cross_modes().is_empty() {
        eprintln!("no cross-cell kernel on this CPU: the cross-cell batches are not compared");
        return;
    }
    let mut seed = 0xc055_u64;
    let size = 4096usize;
    let mask = (size - 1) as i32;
    let mut cases = CrossCases::default();
    let mut calls = 0;
    let mut sizes = std::collections::BTreeSet::new();
    for (window, word_size) in [(40, 3), (16, 2), (40, 1), (3, 5)] {
        for diag_offset in DIAG_OFFSETS {
            for (round, cell_count) in [7usize, 40, 90, 200, 13, 120].into_iter().enumerate() {
                let test = XDiagTest {
                    diag_offset,
                    window,
                    word_size,
                    prefetch: 0,
                };
                let s0 = (next(&mut seed) % 4000) as i32;
                let entries = entries_near(&mut seed, size, s0, test);
                let mut s = s0;
                let cells: Vec<(i32, Vec<i32>)> = (0..cell_count)
                    .map(|k| {
                        let n = match next(&mut seed) % 8 {
                            0 => 25 + next(&mut seed) as usize % 96,
                            1 if round % 3 == 2 && k == cell_count / 2 => 333,
                            _ => next(&mut seed) as usize % 41,
                        };
                        sizes.insert(n);
                        let base = (next(&mut seed) % 3) as i32 * 7;
                        let cell = (s, distinct_offsets(&mut seed, size, n, base));
                        s += 1 + (next(&mut seed) % 4 / 3 * (next(&mut seed) % 6)) as i32;
                        cell
                    })
                    .collect();
                let reference = reference_calls(&entries, mask, test, &cells);
                for chunk in [X_CROSS_TEST_CHUNK, X_CROSS_CHUNK] {
                    let windows = cases.windows;
                    cases.add(&cells, mask, &reference, chunk);
                    assert_eq!(
                        cross_splits(&entries, mask, test, &cells, chunk),
                        (cases.windows - windows) as u64,
                        "chunk {chunk} {test:?}: windows with a repeated entry"
                    );
                }
                calls += same_in_every_mode(&entries, mask, test, &cells);
            }
        }
    }
    assert!(calls > 1000, "the cells must give candidates: {calls}");
    assert!((0..=40).all(|n| sizes.contains(&n)), "{sizes:?}");
    assert!(cases.windows > 0, "{cases:?}");
    assert!(cases.candidate_before_later_cell > 100, "{cases:?}");
    assert!(
        cases.long_cells > 100 && cases.longer_cells > 20,
        "{cases:?}"
    );
    assert!(cases.drains > 50 && cases.drains_in_a_cell > 0, "{cases:?}");
    assert!(cases.moved.iter().all(|&n| n > 0), "{cases:?}");
    assert!(cases.left.iter().all(|&n| n > 0), "{cases:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
// Cells of 1 to 5 hits at consecutive subject offsets, where about half the hits of a cell
// continue a hit of the cell before on its diagonal (q + 1 at s + 1): eight staged hits then
// often hold hits of several cells on one entry, and the window is tested in parts, cut before
// each lane that repeats an entry of its part (parts of every length 1-7, parts that start past
// the window's first lane; with word size 1 the second hit of a diagonal is a candidate, and its
// call's write is what the third reads). Both chunks run (64 drains every few cells).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[test]
fn cross_batches_with_hits_of_one_diagonal_in_consecutive_cells() {
    if cross_modes().is_empty() {
        eprintln!("no cross-cell kernel on this CPU: the cross-cell batches are not compared");
        return;
    }
    let mut seed = 0xd1a9_u64;
    let size = 2048usize;
    let mask = (size - 1) as i32;
    let mut cases = CrossCases::default();
    let mut calls = 0;
    for (window, word_size) in [(40, 1), (40, 3), (4, 1), (2, 1)] {
        for diag_offset in DIAG_OFFSETS {
            for _ in 0..6 {
                let test = XDiagTest {
                    diag_offset,
                    window,
                    word_size,
                    prefetch: 0,
                };
                let s0 = (next(&mut seed) % 4000) as i32;
                let entries = entries_near(&mut seed, size, s0, test);
                let mut cells: Vec<(i32, Vec<i32>)> = Vec::new();
                for k in 0..300 {
                    let n = 1 + next(&mut seed) as usize % 5;
                    let mut hits: Vec<i32> = Vec::new();
                    if let Some((_, before)) = cells.last() {
                        for &q in before {
                            if hits.len() < n && next(&mut seed).is_multiple_of(2) {
                                // The next offset of the diagonal, at most the mask (as the
                                // lookup's offsets are).
                                hits.push((q + 1) & mask);
                            }
                        }
                    }
                    while hits.len() < n {
                        let q = (next(&mut seed) % size as u32) as i32;
                        if !hits.contains(&q) {
                            hits.push(q);
                        }
                    }
                    // The cell's order: the continued hits are not always first.
                    let rotate = next(&mut seed) as usize % n;
                    hits.rotate_left(rotate);
                    cells.push((s0 + k, hits));
                }
                let reference = reference_calls(&entries, mask, test, &cells);
                for chunk in [X_CROSS_TEST_CHUNK, X_CROSS_CHUNK] {
                    let windows = cases.windows;
                    cases.add(&cells, mask, &reference, chunk);
                    assert_eq!(
                        cross_splits(&entries, mask, test, &cells, chunk),
                        (cases.windows - windows) as u64,
                        "chunk {chunk} {test:?}: windows with a repeated entry"
                    );
                }
                calls += same_in_every_mode(&entries, mask, test, &cells);
            }
        }
    }
    assert!(calls > 1000, "the cells must give candidates: {calls}");
    assert!(cases.windows > 1000, "{cases:?}");
    assert!(cases.cuts[1..].iter().all(|&n| n > 0), "{cases:?}");
    assert!(cases.late_parts > 100, "{cases:?}");
    assert!(cases.candidate_before_later_cell > 100, "{cases:?}");
    assert!(cases.drains > 100, "{cases:?}");
}

// No NCBI counterpart: LOSAT's stage; the checks of its precondition have no NCBI form.
/// Runs [`x_scan_cells_cross`] (with `w128`, LOSAT_X_BXCROSS128's [`x_scan_cells_cross128`]) on
/// a one-protein lookup whose `query_length` is set to `query_length` over a table of `mask +
/// 1` entries (the scan stops at its first assert).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn scan_cross_with_bounds(w128: bool, query_length: i32, entries: usize, mask: i32) {
    use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
    use crate::algorithm::tblastx::translation::QueryFrame;
    use crate::config::ScoringMatrix;
    use crate::stats::lookup_protein_params_ungapped;

    let seq = [1u8, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17];
    let mut aa_seq = vec![0u8];
    aa_seq.extend_from_slice(&seq);
    aa_seq.push(0);
    let frames = vec![vec![QueryFrame {
        frame: 1,
        aa_seq,
        aa_seq_nomask: None,
        aa_len: seq.len(),
        orig_len: seq.len(),
        seg_masks: Vec::new(),
    }]];
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (mut lookup, _) = build_ncbi_lookup(&frames, 11, &params, false);
    lookup.query_length = query_length;
    let mut entries = vec![DiagStruct::clear(40); entries];
    let mut table = XDiagTable::new(&mut entries, mask);
    let test = XDiagTest {
        diag_offset: 40,
        window: 40,
        word_size: 3,
        prefetch: 0,
    };
    let subject = [1u8, 3, 4, 5, 6, 7, 0];
    let mut candidate = |_: i32, _: i32, _: i32, _: i32, _: &mut DiagStruct| {};
    #[cfg(target_arch = "x86_64")]
    if !w128 {
        x_scan_cells_cross::<X_CROSS_TEST_CHUNK, _, _>(
            false,
            false,
            &lookup,
            &subject,
            subject.len() - 1,
            &mut table,
            test,
            &mut candidate,
        );
        return;
    }
    assert!(w128, "the AVX2 scan runs on x86_64 only");
    x_scan_cells_cross128::<X_CROSS_TEST_CHUNK, _, _>(
        false,
        &lookup,
        &subject,
        subject.len() - 1,
        &mut table,
        test,
        &mut candidate,
    );
}

// The scan asserts, once per subject and before it stages a hit, that every query offset of
// the lookup is at most the mask: here the lookup is longer than the table.
#[cfg(target_arch = "x86_64")]
#[test]
#[should_panic(expected = "needs every query offset at most the diagonal mask")]
fn cross_scan_asserts_that_the_lookup_fits_the_table() {
    if !x_cross_available() {
        // The scan would panic on the CPU check with another message: fail as expected.
        panic!("needs every query offset at most the diagonal mask (no AVX2 and POPCNT: skipped)");
    }
    scan_cross_with_bounds(false, 64, 16, 15);
}

// No NCBI counterpart: LOSAT's stage; the checks of its precondition have no NCBI form.
// LOSAT_X_BXCROSS128: the 128-bit scan makes the same check (the check is
// `x_cross_check_bounds`, the AVX2 scan's; libtest ignores a `should_panic` test on wasm32,
// whose panics abort, so this runs natively only).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[test]
#[should_panic(expected = "needs every query offset at most the diagonal mask")]
fn cross128_scan_asserts_that_the_lookup_fits_the_table() {
    if !x_cross128_available() {
        panic!(
            "needs every query offset at most the diagonal mask (no SSSE3, SSE4.1 and POPCNT: skipped)"
        );
    }
    scan_cross_with_bounds(true, 64, 16, 15);
}

// The scan also asserts that `mask + 1` is a power of two: the recovery of the query offset
// needs it (here the table has 21 entries).
#[cfg(target_arch = "x86_64")]
#[test]
#[should_panic(expected = "needs a diagonal mask of the form 2^k - 1")]
fn cross_scan_asserts_that_the_mask_is_a_power_of_two_less_one() {
    if !x_cross_available() {
        panic!("needs a diagonal mask of the form 2^k - 1 (no AVX2 and POPCNT: skipped)");
    }
    scan_cross_with_bounds(false, 5, 21, 20);
}

// No NCBI counterpart: LOSAT's stage; the checks of its precondition have no NCBI form.
// LOSAT_X_BXCROSS128: the 128-bit scan makes the same check (natively only, as above).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[test]
#[should_panic(expected = "needs a diagonal mask of the form 2^k - 1")]
fn cross128_scan_asserts_that_the_mask_is_a_power_of_two_less_one() {
    if !x_cross128_available() {
        panic!("needs a diagonal mask of the form 2^k - 1 (no SSSE3, SSE4.1 and POPCNT: skipped)");
    }
    scan_cross_with_bounds(true, 5, 21, 20);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
// LOSAT_X_BXCROSSDEFER's two rare paths, on tables made for them, against the per-hit test
// (every mode, both chunks, [`same_in_every_mode`]), counted by the test probe of the deferred
// arm: cells at subject offsets 5 apart on the same diagonals (q at s, q + 5 at s + 5), each
// diagonal's entry first 5 behind, so a cell's hits are candidates and the next cell's batch
// (cells of eight hits) or the next part of the window (cells of four: two cells per window,
// on the same entries) meets the marker of their deferred calls; and entries that hold the
// marker's value (flag 1, last_hit -1, as a call near offsets of 2^31 stores) under a drain's
// first batch or window, nothing deferred, and under a batch or a window's first part after a
// cell of candidates (the calls made, then the unit again, nothing deferred).
#[cfg(target_arch = "x86_64")]
#[test]
fn deferred_calls_meet_their_marker_and_a_stored_marker_is_a_value() {
    // Only the AVX2 kernel defers its calls (LOSAT_X_BXCROSS128's does not).
    if !x_cross_available() {
        eprintln!("no AVX2 and POPCNT on this CPU: the deferred calls are not compared");
        return;
    }
    let size = 1024usize;
    let mask = (size - 1) as i32;
    let test = XDiagTest {
        diag_offset: 40,
        window: 40,
        word_size: 3,
        prefetch: 0,
    };
    let s0 = 100;
    // `count` cells of `hits` hits on the diagonals 7 + 3 j (j below `hits`) from `s0`,
    // subject offsets `step` apart, and a table whose diagonals last hit 5 before `s0`.
    let cells_on = |count: i32, hits: i32, step: i32, s0: i32| -> Vec<(i32, Vec<i32>)> {
        (0..count)
            .map(|k| {
                let s = s0 + step * k;
                (s, (0..hits).map(|j| (7 + 3 * j + s) & mask).collect())
            })
            .collect()
    };
    let behind = |hits: i32, s0: i32| {
        let mut entries = vec![DiagStruct::default(); size];
        for j in 0..hits {
            let last = s0 - 5 + test.diag_offset;
            entries[(7 + 3 * j) as usize] = entry_of(last as u32 & X_LAST_HIT_MASK);
        }
        entries
    };
    let marker = entry_of(X_CROSS_MARKER as u32);
    // The candidate calls and the probe's [met with calls deferred, met with none].
    let run = |entries: &[DiagStruct], cells: &[(i32, Vec<i32>)]| {
        X_CROSS_MARKS_MET.set([0; 2]);
        let calls = same_in_every_mode(entries, mask, test, cells);
        (calls, X_CROSS_MARKS_MET.get())
    };
    for hits in [8, 4] {
        let (calls, met) = run(&behind(hits, s0), &cells_on(40, hits, 5, s0));
        assert!(calls > 40, "cells of {hits}: {calls} calls");
        assert!(met[0] > 10, "cells of {hits}: {met:?}");
        assert_eq!(met[1], 0, "cells of {hits}: {met:?}");
    }
    // A stored marker under the first batch or window of a drain, then under a batch after a
    // cell of candidates (on other diagonals: 8 + 3 j).
    for hits in [8, 4] {
        let mut entries = vec![DiagStruct::default(); size];
        entries[7] = marker;
        let (calls, met) = run(&entries, &cells_on(2, hits, 5, s0));
        assert!(
            met[1] >= 2 && met[0] == 0,
            "cells of {hits}: {calls} calls, {met:?}"
        );
    }
    let mut entries = behind(8, s0);
    entries[8] = marker;
    let mut cells = cells_on(1, 8, 5, s0);
    cells.push((
        s0 + 1,
        (0..8).map(|j| (8 + 3 * j + s0 + 1) & mask).collect(),
    ));
    let (calls, met) = run(&entries, &cells);
    assert!(calls >= 8, "{calls} calls");
    assert!(met[0] >= 2 && met[1] >= 2, "{met:?}");
    // The same with two cells of four on the diagonals 100 + 3 j after the cell of candidates:
    // one window, whose first part meets the stored marker at entry 100.
    let mut entries = behind(8, s0);
    entries[100] = marker;
    let mut cells = cells_on(1, 8, 5, s0);
    for s in [s0 + 1, s0 + 6] {
        cells.push((s, (0..4).map(|j| (100 + 3 * j + s) & mask).collect()));
    }
    let (calls, met) = run(&entries, &cells);
    assert!(calls >= 8, "{calls} calls");
    assert!(met[0] >= 2 && met[1] >= 2, "{met:?}");
}

// Every pattern of equal lanes of eight coordinates (the 4,140 partitions of the lanes, as
// restricted growth strings): the conflict test sees whether two lanes are equal, and the part
// from each lane ends at the first later lane whose coordinate a lane of the part holds.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[test]
fn the_conflict_test_and_the_parts_follow_the_equal_lanes() {
    if cross_modes().is_empty() {
        eprintln!("no cross-cell kernel on this CPU: the conflict test is not run");
        return;
    }
    // Every kernel's conflict test this machine runs (LOSAT_X_BXCROSS's AVX2 kernel,
    // LOSAT_X_BXCROSS128's) must give the same answer.
    let conflict = |coords: [i32; X_DIAG_LANES]| {
        let mut answers = Vec::new();
        #[cfg(target_arch = "x86_64")]
        if x_cross_available() {
            // SAFETY: the CPU has AVX2 and POPCNT; the load reads the eight lanes of `coords`.
            answers.push(unsafe { XCrossAvx2::conflict(XCrossAvx2::load8(coords.as_ptr())) });
        }
        if x_cross128_available() {
            // SAFETY: the CPU has the 128-bit kernel's instructions; the loads read the eight
            // lanes of `coords`.
            answers.push(unsafe { XCrossW128::conflict(XCrossW128::load8(coords.as_ptr())) });
        }
        assert!(!answers.is_empty(), "no cross-cell kernel on this CPU");
        assert!(
            answers.iter().all(|&a| a == answers[0]),
            "{coords:?}: {answers:?}"
        );
        answers[0]
    };
    let mut classes = [0i32; X_DIAG_LANES];
    let mut patterns = 0;
    loop {
        let coords = classes.map(|class| class * 4_099 + 17);
        let repeated = (1..X_DIAG_LANES).any(|j| coords[..j].contains(&coords[j]));
        assert_eq!(conflict(coords), repeated, "{classes:?}");
        // The parts are the AVX2 kernel's (the 128-bit kernel tests such a window hit by hit).
        #[cfg(target_arch = "x86_64")]
        for start in 0..X_DIAG_LANES {
            // The part's classes so far; it ends at the first class it already holds.
            let mut held = [false; X_DIAG_LANES];
            let mut end = start;
            while end < X_DIAG_LANES && !held[classes[end] as usize] {
                held[classes[end] as usize] = true;
                end += 1;
            }
            assert_eq!(x_cross_part_end(&coords, start), end, "{classes:?} {start}");
        }
        patterns += 1;
        // The next restricted growth string: lane i takes a class up to one above the
        // largest of the lanes before it.
        let Some(i) = (1..X_DIAG_LANES)
            .rev()
            .find(|&i| classes[i] <= *classes[..i].iter().max().unwrap())
        else {
            break;
        };
        classes[i] += 1;
        classes[i + 1..].fill(0);
    }
    assert_eq!(patterns, 4_140);
    assert!(!conflict([3, 1, 2, 0, 7, 6, 5, 4096]));
    assert!(conflict([9; X_DIAG_LANES]));
}

// No NCBI counterpart: LOSAT's choice of its cell test.
// EXPERIMENT (LOSAT_X_BXDIAG, LOSAT_X_BXCROSS, LOSAT_X_BXCROSS128): the cell test of blastx's
// fused loop for each class of CPU and each value of the switches, on whatever machine runs the
// test (the integration tests see only the class of the machine they run on): x86-64 with AVX2
// and POPCNT, with AVX2 but no POPCNT, with SSSE3, SSE4.1 and POPCNT but no AVX2, with none of
// them, and wasm32 simd128; a cell test asked for by name on each.
#[test]
fn the_cell_test_follows_the_cpu_and_the_switches() {
    let per_hit = XCellTest::Diag(XDiagMode::PerHit);
    let none = XCellCpu {
        avx2: false,
        cross: false,
        cross128: false,
    };
    // (CPU, LOSAT_X_BXCROSS on, LOSAT_X_BXCROSSDEFER on, LOSAT_X_BXCROSS128's value, the test)
    #[allow(unused_mut)]
    let mut table = vec![
        (none, true, true, 1, per_hit),
        (none, true, true, 2, per_hit),
    ];
    #[cfg(target_arch = "x86_64")]
    {
        let masked = XCellTest::Diag(XDiagMode::Avx2DeferMasked);
        let avx2 = XCellCpu {
            avx2: true,
            cross: true,
            cross128: true,
        };
        // POPCNT is needed by both kernels.
        let avx2_no_popcnt = XCellCpu {
            avx2: true,
            cross: false,
            cross128: false,
        };
        let sse41 = XCellCpu {
            avx2: false,
            cross: false,
            cross128: true,
        };
        table.extend([
            (avx2, true, true, 1, XCellTest::Cross { defer: true }),
            (avx2, true, false, 1, XCellTest::Cross { defer: false }),
            (avx2, true, true, 0, XCellTest::Cross { defer: true }),
            (avx2, true, true, 2, XCellTest::Cross128),
            (avx2, true, false, 2, XCellTest::Cross128),
            (avx2, false, true, 1, masked),
            (avx2, false, true, 2, masked),
            (avx2_no_popcnt, true, true, 1, masked),
            (avx2_no_popcnt, true, true, 2, masked),
            (sse41, true, true, 1, XCellTest::Cross128),
            (sse41, true, false, 1, XCellTest::Cross128),
            (sse41, true, true, 2, XCellTest::Cross128),
            (sse41, true, true, 0, per_hit),
            (sse41, false, true, 1, per_hit),
            (sse41, false, true, 2, per_hit),
        ]);
    }
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    {
        let simd128 = XCellCpu {
            avx2: false,
            cross: false,
            cross128: true,
        };
        let this = XCellCpu::this();
        assert!(!this.avx2 && !this.cross && this.cross128, "{this:?}");
        table.extend([
            (simd128, true, true, 1, XCellTest::Cross128),
            (simd128, true, false, 1, XCellTest::Cross128),
            (simd128, true, true, 2, XCellTest::Cross128),
            (simd128, true, true, 0, per_hit),
            (simd128, false, true, 1, per_hit),
        ]);
    }
    for &(cpu, cross, defer, cross128, expected) in &table {
        assert_eq!(
            x_cell_test_on("", cross, defer, cross128, cpu),
            expected,
            "{cpu:?} cross {cross} defer {defer} cross128 {cross128}"
        );
        for (name, mode) in [
            ("per-hit", XDiagMode::PerHit),
            ("scalar", XDiagMode::Scalar),
        ] {
            assert_eq!(
                x_cell_test_on(name, cross, defer, cross128, cpu),
                XCellTest::Diag(mode),
                "{name} {cpu:?}"
            );
        }
    }
}

// EXPERIMENT (LOSAT_X_BXCROSS): the query offset recovered from the entry coordinate and the
// subject offset is the hit's for every query offset up to the mask, with subject offsets
// below, at and past it, up to i32::MAX (the sum wraps).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[test]
fn the_query_offset_is_recovered_from_the_coordinate() {
    for mask in [0i32, 1, 7, 1_023, 65_535, (1 << 30) - 1] {
        let queries = [0, 1, mask / 3, mask / 2, mask - 1, mask];
        let past = mask.saturating_mul(2).saturating_add(3);
        let subjects = [0, 1, 5, mask, mask + 1, past, 5_000, i32::MAX - 1, i32::MAX];
        for q in queries.into_iter().filter(|q| (0..=mask).contains(q)) {
            for s in subjects {
                let coord = q.wrapping_sub(s) & mask;
                assert_eq!(x_cross_query(coord, s, mask), q, "q {q} s {s} mask {mask}");
            }
        }
    }
}
