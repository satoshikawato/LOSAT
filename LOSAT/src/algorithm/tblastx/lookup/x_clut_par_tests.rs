use super::super::tests::{queries, XorShift, FEW_LETTERS, LETTERS};
use super::*;

// Item SC9i of line NP: a one-latch pass makes every item once, in order, at any thread count,
// and an item that panics on a spawned task or on the calling thread ends the pass with its
// panic instead of leaving the caller waiting.
#[test]
fn one_latch_pass_makes_each_item_once_and_passes_a_panic_on() {
    for threads in [1usize, 2, 3, 8] {
        for (items, extras) in [(0usize, 0usize), (1, 0), (5, 3), (200, 17)] {
            let extra_runs: Vec<AtomicUsize> = (0..extras).map(|_| AtomicUsize::new(0)).collect();
            let made = in_pool(threads, || {
                x_one_latch_pass(
                    threads,
                    items,
                    extras,
                    || 0usize,
                    |taken, index| {
                        *taken += 1;
                        index * 3
                    },
                    |index| {
                        extra_runs[index].fetch_add(1, Ordering::Relaxed);
                    },
                )
            });
            assert_eq!(made, (0..items).map(|index| index * 3).collect::<Vec<_>>());
            assert!(extra_runs
                .iter()
                .all(|runs| runs.load(Ordering::Relaxed) == 1));
        }
        for bad in [0usize, 7, 63] {
            let outcome = std::panic::catch_unwind(|| {
                in_pool(threads, || {
                    x_one_latch_pass(
                        threads,
                        64,
                        0,
                        || (),
                        |_, index| {
                            assert_ne!(index, bad, "item {bad} fails");
                            index
                        },
                        |_| {},
                    )
                })
            });
            assert!(
                outcome.is_err(),
                "threads {threads}: item {bad}'s panic passed on"
            );
        }
    }
}

// Item SC9i of line NP: a gang's passes make every item of each pass once, in order, at any
// thread count, one pass after another on the same helpers; an item that panics on a helper or
// on the calling thread ends the build with its panic instead of leaving a thread waiting.
#[test]
fn gang_passes_make_each_item_once_and_pass_a_panic_on() {
    for threads in [1usize, 2, 3, 8] {
        let made = in_pool(threads, || {
            x_with_gang(threads, |gang| {
                [0usize, 1, 5, 200, 0, 64, 3]
                    .iter()
                    .enumerate()
                    .map(|(pass, &items)| {
                        let runs: Vec<AtomicUsize> =
                            (0..items).map(|_| AtomicUsize::new(0)).collect();
                        let results = gang.pass(
                            items,
                            pass == 6,
                            || 0usize,
                            |taken, index| {
                                *taken += 1;
                                runs[index].fetch_add(1, Ordering::Relaxed);
                                index * 3 + pass
                            },
                        );
                        assert!(runs.iter().all(|runs| runs.load(Ordering::Relaxed) == 1));
                        results
                    })
                    .collect::<Vec<_>>()
            })
        });
        for (pass, (results, items)) in made.iter().zip([0usize, 1, 5, 200, 0, 64, 3]).enumerate() {
            assert_eq!(
                *results,
                (0..items).map(|index| index * 3 + pass).collect::<Vec<_>>(),
                "threads {threads}, pass {pass}"
            );
        }
        for bad in [0usize, 7, 63] {
            let outcome = std::panic::catch_unwind(|| {
                in_pool(threads, || {
                    x_with_gang(threads, |gang| {
                        let first = gang.pass(64, false, || (), |_, index| index);
                        let second = gang.pass(
                            64,
                            bad != 7,
                            || (),
                            |_, index| {
                                assert_ne!(index, bad, "item {bad} fails");
                                index
                            },
                        );
                        (first, second)
                    })
                })
            });
            assert!(
                outcome.is_err(),
                "threads {threads}: item {bad}'s panic passed on"
            );
        }
    }
}

// Item SC9i of line NP: helpers that start after a gang's last pass is posted join it while it
// is open (a build of one wave has only that pass): the caller's first item waits until another
// thread has made one (a minute at most).
#[test]
fn helpers_join_a_gangs_last_pass_while_it_is_open() {
    for threads in [2usize, 4] {
        let others = AtomicUsize::new(0);
        let joined = in_pool(threads, || {
            let caller = std::thread::current().id();
            x_with_gang(threads, |gang| {
                let made = gang.pass(
                    64,
                    true,
                    || (),
                    |_, index| {
                        if std::thread::current().id() != caller {
                            others.fetch_add(1, Ordering::AcqRel);
                        } else if index == 0 {
                            let deadline =
                                std::time::Instant::now() + std::time::Duration::from_secs(60);
                            while others.load(Ordering::Acquire) == 0
                                && std::time::Instant::now() < deadline
                            {
                                std::thread::yield_now();
                            }
                        }
                        index
                    },
                );
                assert_eq!(made, (0..64).collect::<Vec<_>>());
                others.load(Ordering::Acquire) > 0
            })
        });
        assert!(
            joined,
            "threads {threads}: no helper made an item of the last pass"
        );
    }
}

// Item SC9i of line NP: the overflow cells read from a cell's count are those its appends take
// one by one (the serial count pass's statement).
#[test]
fn overflow_cells_of_a_count_are_those_its_appends_take() {
    let mut taken = 0usize;
    for n in 0..10_000u32 {
        assert_eq!(x_overflow_cells_of(n), taken, "n {n}");
        taken += usize::from(x_takes_overflow_cell(n));
    }
    let counts = [0u32, 1, 5, 6, 7, 8, 11, 12, 40];
    assert_eq!(x_overflow_cells_in(&counts), 1 + 1 + 2 + 2 + 3 + 10);
}

fn in_pool<R: Send>(threads: usize, work: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("test pool")
        .install(work)
}

#[allow(clippy::too_many_arguments)]
fn parallel(
    threads: usize,
    threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    capacity: usize,
    max_appends: Option<usize>,
    chunk_sizes: &[usize],
    stripe_shift: Option<usize>,
    passes: XPasses<'_>,
    stripe_fill: bool,
    spare: Option<(&XFlatSpare, bool)>,
) -> XFlatBuild {
    if let Some((spare, true)) = spare {
        poison(spare);
    }
    in_pool(threads, || {
        x_try_build_flat_parallel_with(
            5,
            threshold,
            query,
            locations,
            tables,
            XParallelCuts {
                overflow_cell_capacity: capacity,
                max_appends,
                threads,
                chunk_sizes,
                stripe_shift,
                // Item SC9i of line NP: one wave, or many (word size 5's 759,376 cells in
                // waves of 2^12 to 2^16 cells; at most one per stripe).
                wave_cells: [None, Some(1 << 12), Some(5003), Some(1 << 16)]
                    [(threads + chunk_sizes.len() + threshold as usize) % 4],
                appends_segment: X_APPENDS_SEGMENT,
                piece_appends: X_PIECE_APPENDS,
                segments: None,
                passes,
                spare: spare.map(|(spare, _)| spare),
                lent: None,
                x_ahead: false,
                stripe_fill,
            },
        )
    })
    .expect("word size 5")
}

/// Fills `spare`'s vectors, to their whole capacity, with `u32::MAX` / `i32::MIN`: an entry
/// that no stripe writes then keeps a value that no build makes (the build takes the vectors'
/// values as they are, module header).
fn poison(spare: &XFlatSpare) {
    let mut kept = spare.0.lock().unwrap();
    let (starts, hits) = &mut *kept;
    starts.resize(starts.capacity(), u32::MAX);
    starts.fill(u32::MAX);
    hits.resize(hits.capacity(), i32::MIN);
    hits.fill(i32::MIN);
}

// The pass implementations (LOSAT_X_BPLKPFLAT): `par_iter`, and one-latch passes with the stripes
// of `x_stripe_build` or (LOSAT_X_BPSTRIPEFILL) of `x_stripe_fill`.
const FORMS: [(XPasses<'static>, bool); 3] = [
    (XPasses::Rayon, false),
    (XPasses::OneLatch(None), false),
    (XPasses::OneLatch(None), true),
];

/// Extra tasks that count how often each one ran.
struct Counted(Vec<AtomicUsize>);

impl Counted {
    fn new(count: usize) -> Self {
        Self((0..count).map(|_| AtomicUsize::new(0)).collect())
    }

    fn assert_each_ran(&self, times: usize, what: &str) {
        for (index, runs) in self.0.iter().enumerate() {
            assert_eq!(
                runs.load(Ordering::Relaxed),
                times,
                "{what}: extra task {index}"
            );
        }
    }
}

impl XExtraTasks for Counted {
    fn count(&self) -> usize {
        self.0.len()
    }

    fn run(&self, index: usize) {
        self.0[index].fetch_add(1, Ordering::Relaxed);
    }
}

/// A spare of garbage: `starts` with room for a word-size-5 table (759,376 cells and one), `hits`
/// of random length.
fn garbage_spare(rng: &mut XorShift) -> XFlatSpare {
    let spare = XFlatSpare::default();
    *spare.0.lock().unwrap() = (
        vec![u32::MAX - 1; 759_377 + rng.below(9)],
        vec![i32::MIN + 1; rng.below(4000)],
    );
    spare
}

/// Reserves room in a spare of garbage for a word-size-5 table's `starts` and for 65,536 hits,
/// so that `poison` fills whole vectors and the build never allocates a fresh (zeroed) vector
/// where a spare is poisoned.
fn give_room(spare: &XFlatSpare) {
    let mut kept = spare.0.lock().unwrap();
    let (starts, hits) = &mut *kept;
    starts.reserve(759_377);
    hits.reserve(1 << 16);
}

/// The capacities of the spare's vectors.
fn spare_capacities(spare: &XFlatSpare) -> (usize, usize) {
    let kept = spare.0.lock().unwrap();
    (kept.0.capacity(), kept.1.capacity())
}

// The cuts of the tests: thread counts, chunk sizes (even and uneven), stripe widths (the
// narrowest, 2^5 cells, a middle one and the widest).
const THREADS: [usize; 4] = [1, 2, 3, 7];
const CHUNK_SIZES: [&[usize]; 4] = [&[], &[1], &[13, 1, 97, 2, 5, 400], &[usize::MAX]];
const SHIFTS: [Option<usize>; 4] = [None, Some(5), Some(11), Some(16)];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
// Random queries (few letters and low thresholds, so that cells overflow; all letters, with
// letters outside the compressed alphabet), SEG-masked runs (several segments per query), several
// NCBI tables and split batches (segments outside every table), threshold 0: at 1, 2, 3 and 7
// threads, with even and uneven chunks and narrow to wide stripes, the parallel table is the
// serial one slot for slot (and the reference's, cell by cell in scan order).
#[test]
fn parallel_table_equals_serial_on_random_queries() {
    let mut rng = XorShift(0x510e_527f_ade6_82d1);
    let mut overflow_seen = false;
    for round in 0..10 {
        let (letters, threshold): (&[u8], i32) = match round % 5 {
            0 => (&FEW_LETTERS[..], 1100),
            1 => (&FEW_LETTERS[..], 1600),
            2 => (&LETTERS[..], 1100),
            3 => (&[1u8, 11, 17, 21, 24, 10, 3][..], 2000),
            _ => (&FEW_LETTERS[..], 0),
        };
        let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 1 + rng.below(300)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 3);
        let mut tables = Vec::new();
        let mut q = 0;
        while q < lengths.len() {
            let end = (q + 1 + rng.below(2)).min(lengths.len());
            if rng.below(3) != 0 {
                tables.push(first[q]..first[end]);
            }
            q = end;
        }
        let serial = flat(x_try_build_flat_compressed_lookup(
            5, threshold, &query, &locations, &tables,
        ));
        let reference =
            try_build_blosum62_compressed_lookup(5, threshold, &query, &locations, &tables)
                .expect("word size 5")
                .expect("not refused");
        x_assert_same_compressed_table(&serial, &reference);
        overflow_seen |= serial.longest_chain > 5;
        // Every thread count with every chunk size and every form (`FORMS`); the stripe
        // widths in turn.
        for (i, threads) in THREADS.into_iter().enumerate() {
            for (j, chunk_sizes) in CHUNK_SIZES.into_iter().enumerate() {
                let stripe_shift = SHIFTS[(i + j + round) % SHIFTS.len()];
                for (passes, stripe_fill) in FORMS {
                    let built = flat(Some(parallel(
                        threads,
                        threshold,
                        &query,
                        &locations,
                        &tables,
                        NCBI_OVERFLOW_CELL_CAPACITY,
                        None,
                        chunk_sizes,
                        stripe_shift,
                        passes,
                        stripe_fill,
                        None,
                    )));
                    x_assert_same_flat_tables(&built, &serial);
                }
            }
        }
    }
    assert!(overflow_seen, "cells with overflow lists");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1310
// Every search matrix but BLOSUM62 and IDENTITY (item SC2b): random queries of few letters, NCBI's
// tables of several batches and split batches, word size 5 with capacities around each table's
// overflow cells, at 2, 3 and 7 threads with even and uneven chunks. The pool's build of the
// matrix's words gives the serial builder's table of those words (`x_try_build_flat_of`) slot
// for slot, which is the reference's cell by cell, or refuses where both refuse; a table of the
// matrix's words is not BLOSUM62's. BLOSUM45 and PAM30 also at word size 6.
#[test]
fn parallel_tables_of_each_matrix_equal_its_serial_tables() {
    let mut rng = XorShift(0x1f83_d9ab_fb41_bd6b);
    let (mut built_seen, mut refused_seen) = (0usize, 0usize);
    let (mut overflow_seen, mut differs_seen) = (false, false);
    let cuts = |threads, capacity, chunk_sizes| XParallelCuts {
        overflow_cell_capacity: capacity,
        max_appends: None,
        threads,
        chunk_sizes,
        stripe_shift: None,
        wave_cells: None,
        appends_segment: X_APPENDS_SEGMENT,
        piece_appends: X_PIECE_APPENDS,
        segments: None,
        passes: XPasses::Rayon,
        spare: None,
        lent: None,
        x_ahead: false,
        stripe_fill: false,
    };
    for matrix in ScoringMatrix::ALL {
        if matches!(matrix, ScoringMatrix::Blosum62 | ScoringMatrix::Identity) {
            continue;
        }
        let compressed = CompressedScores::score_matrix(matrix, 15).expect("compressed matrix");
        for round in 0..3 {
            let threshold = [1100, 1600, 2000][round];
            let lengths: Vec<usize> = (0..2 + rng.below(4)).map(|_| 5 + rng.below(250)).collect();
            let (query, locations, first) = queries(&mut rng, &lengths, &FEW_LETTERS, round % 2);
            let mut tables = Vec::new();
            let mut q = 0;
            while q < lengths.len() {
                let end = (q + 1 + rng.below(3)).min(lengths.len());
                if rng.below(4) != 0 {
                    tables.push(first[q]..first[end]);
                }
                q = end;
            }
            let mut capacities = vec![1usize, NCBI_OVERFLOW_CELL_CAPACITY];
            for table in &tables {
                let mut own = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
                own.threshold = threshold;
                own.overflow_cell_capacity = usize::MAX;
                own.add_neighboring_words(&compressed, &query, &locations[table.clone()]);
                let cells = own.overflow_cells.len();
                capacities.extend([cells.saturating_sub(1), cells, cells + 1]);
            }
            capacities.retain(|&capacity| capacity > 0);
            capacities.sort_unstable();
            capacities.dedup();
            for capacity in capacities {
                let what = format!("{matrix:?} threshold {threshold} capacity {capacity}");
                let reference = try_build_compressed_lookup_with_capacity(
                    matrix, 5, threshold, &query, &locations, &tables, capacity,
                )
                .expect("word size 5");
                let serial = x_try_build_flat_of(
                    matrix,
                    5,
                    threshold,
                    &query,
                    &locations,
                    &tables,
                    capacity,
                    X_COUNT_CHUNK,
                )
                .expect("word size 5");
                match (&reference, &serial) {
                    (Ok(reference), XFlatBuild::Table(serial)) => {
                        x_assert_same_compressed_table(serial, reference);
                        overflow_seen |= serial.longest_chain > 5;
                        built_seen += 1;
                    }
                    (Err(want), XFlatBuild::Refused(got)) if want == got => refused_seen += 1,
                    _ => panic!("{what}: the serial outcome is not the reference's"),
                }
                for (threads, chunk_sizes) in [(2, &[][..]), (3, &[7, 1, 50][..]), (7, &[1][..])] {
                    let built = in_pool(threads, || {
                        x_try_build_flat_parallel_of(
                            matrix,
                            5,
                            threshold,
                            &query,
                            &locations,
                            &tables,
                            cuts(threads, capacity, chunk_sizes),
                        )
                    });
                    match (&serial, built) {
                        (XFlatBuild::Table(serial), Some(XFlatBuild::Table(built))) => {
                            x_assert_same_flat_tables(&built, serial)
                        }
                        (XFlatBuild::Refused(got), Some(XFlatBuild::Refused(want))) if *got == want => {}
                        _ => panic!("{what} n{threads}: the pool's outcome is not the serial one"),
                    }
                }
            }
            // The matrix's words, not BLOSUM62's, on the same query.
            if let (Some(XFlatBuild::Table(ours)), Some(XFlatBuild::Table(blosum62))) = (
                x_try_build_flat_of(
                    matrix,
                    5,
                    threshold,
                    &query,
                    &locations,
                    &[],
                    NCBI_OVERFLOW_CELL_CAPACITY,
                    X_COUNT_CHUNK,
                ),
                x_try_build_flat_compressed_lookup(5, threshold, &query, &locations, &[]),
            ) {
                differs_seen |= ours.hits != blosum62.hits;
            }
        }
    }
    // Word size 6 (a backbone of 15^6 cells) with two matrices, NCBI's capacity.
    for matrix in [ScoringMatrix::Blosum45, ScoringMatrix::Pam30] {
        let lengths: Vec<usize> = (0..2 + rng.below(3)).map(|_| 1 + rng.below(200)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, &LETTERS[..], 1);
        let tables = vec![first[0]..first[lengths.len()]];
        let Some(XFlatBuild::Table(serial)) = x_try_build_flat_of(
            matrix,
            6,
            2100,
            &query,
            &locations,
            &tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
            X_COUNT_CHUNK,
        ) else {
            panic!("{matrix:?} word size 6: no serial table")
        };
        for (threads, chunk_sizes) in [(2, &[][..]), (3, &[7, 50][..]), (7, &[1][..])] {
            let built = in_pool(threads, || {
                x_try_build_flat_parallel_of(
                    matrix,
                    6,
                    2100,
                    &query,
                    &locations,
                    &tables,
                    cuts(threads, NCBI_OVERFLOW_CELL_CAPACITY, chunk_sizes),
                )
            });
            let Some(XFlatBuild::Table(built)) = built else {
                panic!("{matrix:?} word size 6 n{threads}: no parallel table")
            };
            x_assert_same_flat_tables(&built, &serial);
        }
    }
    assert!(
        built_seen > 0 && refused_seen > 0 && overflow_seen && differs_seen,
        "built {built_seen} refused {refused_seen} overflow {overflow_seen} differs {differs_seen}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_try_build_flat_in_pool`)
// The entry of blastp-fast with another matrix (`x_bp_build_compressed_lookup_of`, items SC2b,
// SC2d-2): the serial builders at one thread and the pool's build at 3 threads give the matrix's
// serial table, with rayon's passes and with the one-latch passes (LOSAT_X_BPLKPFLAT) on a spare
// of garbage, with the stripe fill's switch as it is set; a build ahead of its block
// (LOSAT_X_BPLKPAHEAD) gives it at once and in its two parts (`x_bp_flat_ahead_appends`,
// `x_bp_flat_ahead_table`) on the vectors `x_make_room` made ready, holding them, and gives
// nothing without a spare. The two passes apart refuse with the matrix's words as the build at
// once does (a capacity of one cell).
#[test]
fn the_entry_of_another_matrix_builds_its_serial_table_at_any_thread_count() {
    let mut rng = XorShift(0x6a09_e667_f3bc_c908);
    for matrix in [ScoringMatrix::Blosum45, ScoringMatrix::Pam250] {
        let lengths: Vec<usize> = (0..3).map(|_| 20 + rng.below(200)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, &LETTERS[..], 1);
        let tables = vec![first[0]..first[lengths.len()]];
        let Some(XFlatBuild::Table(serial)) = x_try_build_flat_of(
            matrix,
            5,
            1100,
            &query,
            &locations,
            &tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
            X_COUNT_CHUNK,
        ) else {
            panic!("{matrix:?}: no serial table")
        };
        for threads in [1, 3] {
            let built = crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                Ok(super::super::x_bp_build_compressed_lookup_of(
                    matrix,
                    5,
                    1100,
                    &query,
                    &locations,
                    &tables,
                    pool,
                    XFlatArgs::default(),
                ))
            })
            .expect("pool");
            let Some(Ok(XCompressedLookup::Flat(built))) = built else {
                panic!("{matrix:?} n{threads}: no flat table")
            };
            x_assert_same_flat_tables(&built, &serial);
        }
        let backbone_size = serial.starts.len() - 1;
        let total = serial.hits.len();
        let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
        let build = |args: XFlatArgs<'_>| {
            crate::utils::threading::with_search_pool(3, "blastp", |pool| {
                Ok(super::super::x_bp_build_compressed_lookup_of(
                    matrix, 5, 1100, &query, &locations, &tables, pool, args,
                ))
            })
            .expect("pool")
        };
        // The one-latch passes on a spare of garbage.
        let spare = garbage_spare(&mut rng);
        let Some(Ok(XCompressedLookup::Flat(built))) = build(XFlatArgs {
            mode: XLkpFlat::On,
            spare: Some(&spare),
            extra: None,
            segments: None,
            x_ahead: false,
            appends: 0,
        }) else {
            panic!("{matrix:?}: no flat table on the spare")
        };
        x_assert_same_flat_tables(&built, &serial);
        // A build ahead of its block, at once and in its two parts: every vector from the spare
        // made ready.
        let ready = |rng: &mut XorShift| {
            let spare = XFlatSpare::default();
            let mut hits = Vec::with_capacity(total + 5);
            hits.extend((0..total / 2).map(|_| rng.below(1 << 20) as i32));
            *spare.0.lock().unwrap() = (vec![u32::MAX - 3; backbone_size + 1], hits);
            *spare.1.lock().unwrap() = vec![u32::MAX; pv_words];
            let block = XFlatSpare::default();
            assert_eq!(spare.x_make_room(&block, backbone_size, total, total), 0);
            spare
        };
        let held = |built: &XFlatCompressedLookup| {
            [
                built.starts.as_ptr() as usize,
                built.hits.as_ptr() as usize,
                built.pv.as_ptr() as usize,
            ]
        };
        let spare = ready(&mut rng);
        let given = spare_addresses(&spare);
        let Some(Ok(XCompressedLookup::Flat(built))) = build(XFlatArgs {
            mode: XLkpFlat::On,
            spare: Some(&spare),
            extra: None,
            segments: None,
            x_ahead: true,
            appends: 0,
        }) else {
            panic!("{matrix:?}: no flat table built ahead")
        };
        x_assert_same_flat_tables(&built, &serial);
        assert_eq!(held(&built), given, "{matrix:?}: the spare's vectors");
        let spare = ready(&mut rng);
        let given = spare_addresses(&spare);
        let built = crate::utils::threading::with_search_pool(3, "blastp", |pool| {
            let appends = super::super::x_bp_flat_ahead_appends(
                matrix,
                5,
                1100,
                &query,
                &locations,
                &tables,
                pool,
                XFlatArgs {
                    mode: XLkpFlat::On,
                    spare: None,
                    extra: None,
                    segments: None,
                    x_ahead: true,
                    appends: 0,
                },
            )
            .expect("the appends pass");
            Ok(super::super::x_bp_flat_ahead_table(
                matrix, appends, &query, &locations, pool, &spare, None,
            ))
        })
        .expect("pool")
        .expect("a table built ahead in two parts");
        x_assert_same_flat_tables(&built, &serial);
        assert_eq!(
            held(&built),
            given,
            "{matrix:?}: the spare's vectors, two parts"
        );
        // Without a spare a build ahead makes no table.
        assert!(
            build(XFlatArgs {
                mode: XLkpFlat::On,
                spare: None,
                extra: None,
                segments: None,
                x_ahead: true,
                appends: 0,
            })
            .is_none(),
            "{matrix:?}: a table built ahead without a spare"
        );
        // The two passes apart at a capacity of one cell: the matrix's words make the same
        // outcome as its serial builder (the capacity test makes them again).
        let serial_one = x_try_build_flat_of(
            matrix,
            5,
            1100,
            &query,
            &locations,
            &tables,
            1,
            X_COUNT_CHUNK,
        );
        let cuts = |threads| XParallelCuts::<i32> {
            overflow_cell_capacity: 1,
            max_appends: None,
            threads,
            chunk_sizes: &[],
            stripe_shift: None,
            wave_cells: None,
            appends_segment: X_APPENDS_SEGMENT,
            piece_appends: X_PIECE_APPENDS,
            segments: None,
            passes: XPasses::OneLatch(None),
            spare: None,
            lent: None,
            x_ahead: false,
            stripe_fill: true,
        };
        let stage = in_pool(3, || {
            x_flat_appends_on(
                matrix,
                5,
                1100,
                &query,
                &locations,
                &tables,
                cuts(3),
                |_| None,
            )
        });
        let XFlatStage::Appends(appends) = stage else {
            panic!("{matrix:?}: the appends pass ended the build")
        };
        let apart = in_pool(2, || {
            x_flat_stripes_on(matrix, appends, &query, &locations, None, None, false, None)
        });
        match (&serial_one, &apart) {
            (Some(XFlatBuild::Refused(want)), Some(XFlatBuild::Refused(got))) if want == got => {}
            (Some(XFlatBuild::Table(want)), Some(XFlatBuild::Table(got))) => {
                x_assert_same_flat_tables(got, want)
            }
            _ => panic!("{matrix:?}: capacity 1, the outcomes differ"),
        }
    }
}

fn flat(build: Option<XFlatBuild>) -> XFlatCompressedLookup {
    match build.expect("word size 5") {
        XFlatBuild::Table(table) => table,
        XFlatBuild::Refused(filled) => panic!("refused: {filled:?}"),
        XFlatBuild::TooLarge => panic!("too large"),
    }
}

// A query longer than a chunk's offset span (2^16) in one segment, and in several: the chunks
// are cut where the span ends, and the table is the serial one.
#[test]
fn parallel_table_equals_serial_past_the_offset_span() {
    let mut rng = XorShift(0x9b05_688c_2b3e_6c1f);
    for masked in [0, 3] {
        let (query, locations, first) = queries(&mut rng, &[70_000, 900], &LETTERS, masked);
        let tables = [first[0]..first[2]];
        let serial = flat(x_try_build_flat_compressed_lookup(
            5, 2000, &query, &locations, &tables,
        ));
        for (threads, chunk_sizes) in [(2, &[usize::MAX][..]), (3, &[][..]), (7, &[40_000][..])] {
            for (passes, stripe_fill) in FORMS {
                let built = flat(Some(parallel(
                    threads,
                    2000,
                    &query,
                    &locations,
                    &tables,
                    NCBI_OVERFLOW_CELL_CAPACITY,
                    None,
                    chunk_sizes,
                    None,
                    passes,
                    stripe_fill,
                    None,
                )));
                x_assert_same_flat_tables(&built, &serial);
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
// The refusal at small capacities (around each table's own cells and its cells of the table of
// all batches, as the serial builder's test): at 2, 3 and 7 threads and uneven chunks the
// parallel builder refuses exactly when the reference does, and otherwise builds its table.
// The bound on the appends: `TooLarge` exactly when the appends pass it. Both pass
// implementations; the one-latch passes with extra tasks, each of which runs exactly once, also
// when the build is refused or `TooLarge`, and on a spare of garbage (filled with `poison`'s
// pattern or not), which a refused build consumes and a `TooLarge` build drops; the next build on
// a spare of garbage is still exact.
#[test]
fn parallel_refusal_and_bound_are_the_serial_ones() {
    let mut rng = XorShift(0x1f83_d9ab_fb41_bd6b);
    let mut refused_seen = 0usize;
    let mut built_seen = 0usize;
    for round in 0..8 {
        let threshold = [1100, 1600, 2000][round % 3];
        let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 5 + rng.below(350)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, &FEW_LETTERS, round % 2);
        let mut tables = Vec::new();
        let mut q = 0;
        while q < lengths.len() {
            let end = (q + 1 + rng.below(3)).min(lengths.len());
            if rng.below(4) != 0 {
                tables.push(first[q]..first[end]);
            }
            q = end;
        }
        let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
        let mut capacities = vec![1usize];
        for table in &tables {
            let mut own = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
            own.threshold = threshold;
            own.overflow_cell_capacity = usize::MAX;
            own.add_neighboring_words(&matrix, &query, &locations[table.clone()]);
            let cells = own.overflow_cells.len();
            capacities.extend([cells.saturating_sub(1), cells, cells + 1, 2 * cells]);
        }
        capacities.retain(|&capacity| capacity > 0);
        for capacity in capacities {
            let reference = try_build_compressed_lookup_with_capacity(
                Blosum62Compressed,
                5,
                threshold,
                &query,
                &locations,
                &tables,
                capacity,
            )
            .expect("word size 5");
            for (threads, chunk_sizes) in [(2, &[][..]), (3, &[7, 1, 50][..]), (7, &[1][..])] {
                for (one_latch, stripe_fill) in [(false, false), (true, false), (true, true)] {
                    let extra = Counted::new(1 + rng.below(40));
                    let passes = if one_latch {
                        XPasses::OneLatch(Some(&extra))
                    } else {
                        XPasses::Rayon
                    };
                    let spare = garbage_spare(&mut rng);
                    let poison = rng.below(2) == 0;
                    if poison {
                        give_room(&spare);
                    }
                    match (
                        &reference,
                        parallel(
                            threads,
                            threshold,
                            &query,
                            &locations,
                            &tables,
                            capacity,
                            None,
                            chunk_sizes,
                            None,
                            passes,
                            stripe_fill,
                            one_latch.then_some((&spare, poison)),
                        ),
                    ) {
                        (Err(want), XFlatBuild::Refused(got)) if *want == got => refused_seen += 1,
                        (Ok(reference), XFlatBuild::Table(built)) => {
                            x_assert_same_compressed_table(&built, reference);
                            built_seen += 1;
                        }
                        (reference, _) => panic!(
                            "round {round}, capacity {capacity}, threads {threads}: reference refused {}",
                            reference.is_err()
                        ),
                    }
                    extra.assert_each_ran(
                        usize::from(one_latch),
                        &format!("round {round}, capacity {capacity}, threads {threads}"),
                    );
                    assert_eq!(spare_capacities(&spare) == (0, 0), one_latch);
                }
            }
        }
        // The bound: below, at and above the number of appends (the serial table's hits).
        let serial = flat(x_try_build_flat_compressed_lookup(
            5, threshold, &query, &locations, &tables,
        ));
        let appends = serial.hits.len();
        for (max_appends, too_large) in [
            (appends.saturating_sub(1), appends > 0),
            (appends, false),
            (appends + 1, false),
        ] {
            for (one_latch, stripe_fill) in [(false, false), (true, false), (true, true)] {
                let extra = Counted::new(1 + rng.below(40));
                let passes = if one_latch {
                    XPasses::OneLatch(Some(&extra))
                } else {
                    XPasses::Rayon
                };
                let spare = garbage_spare(&mut rng);
                let poison = rng.below(2) == 0;
                if poison {
                    give_room(&spare);
                }
                let built = parallel(
                    3,
                    threshold,
                    &query,
                    &locations,
                    &tables,
                    NCBI_OVERFLOW_CELL_CAPACITY,
                    Some(max_appends),
                    &[5],
                    None,
                    passes,
                    stripe_fill,
                    one_latch.then_some((&spare, poison)),
                );
                assert_eq!(
                    matches!(built, XFlatBuild::TooLarge),
                    too_large,
                    "round {round}, bound {max_appends}, appends {appends}"
                );
                if let XFlatBuild::Table(built) = built {
                    x_assert_same_flat_tables(&built, &serial);
                }
                extra.assert_each_ran(
                    usize::from(one_latch),
                    &format!("round {round}, bound {max_appends}"),
                );
                // A `TooLarge` build drops the spare's vectors, a table takes them.
                assert_eq!(spare_capacities(&spare) == (0, 0), one_latch);
            }
        }
    }
    assert!(refused_seen > 0 && built_seen > 0);
}

// LOSAT_X_BPLKPFLAT: a build on vectors of garbage (`u32::MAX`, `i32::MIN`, random values) longer
// or shorter than the table's, with and without the room for it, and then on the vectors that
// build gives back, is the fresh build slot for slot, with every form and at 1, 2,
// 3 and 7 threads; where the spare's `hits` has the room, the table holds that buffer.
#[test]
fn builds_on_recycled_vectors_equal_fresh_builds() {
    let mut rng = XorShift(0x6a09_e667_f3bc_c908);
    let spare = XFlatSpare::default();
    let mut reused_seen = 0usize;
    for round in 0..6 {
        let (letters, threshold): (&[u8], i32) = match round % 3 {
            0 => (&FEW_LETTERS[..], 1100),
            1 => (&LETTERS[..], 1100),
            _ => (&FEW_LETTERS[..], 1600),
        };
        let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 1 + rng.below(400)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 3);
        let tables = [first[0]..first[lengths.len()]];
        let serial = flat(x_try_build_flat_compressed_lookup(
            5, threshold, &query, &locations, &tables,
        ));
        let cells = serial.starts.len();
        let total = serial.hits.len();
        for case in 0..3 {
            for threads in THREADS {
                for (passes, stripe_fill) in FORMS {
                    let mut random = |len: usize| -> Vec<u32> {
                        (0..len)
                            .map(|_| rng.below(u32::MAX as usize) as u32)
                            .collect()
                    };
                    let garbage: (Vec<u32>, Vec<i32>) = match case {
                        // Longer than the table's.
                        0 => (vec![u32::MAX; cells + 7], vec![i32::MIN; total + 100]),
                        // Random values; `hits` shorter, with room to grow in place.
                        1 => {
                            let starts = random(cells);
                            let mut hits = Vec::with_capacity(2 * total + 1);
                            hits.extend(random(total / 2).into_iter().map(|value| value as i32));
                            (starts, hits)
                        }
                        // No room: fresh vectors.
                        _ => (Vec::new(), vec![i32::MIN; total.saturating_sub(1)]),
                    };
                    *spare.0.lock().unwrap() = garbage;
                    // The garbage, then the vectors the first build gives back.
                    for build in 0..2 {
                        let (given, room) = {
                            let kept = spare.0.lock().unwrap();
                            (kept.1.as_ptr() as usize, kept.1.capacity() >= total)
                        };
                        let built = flat(Some(parallel(
                            threads,
                            threshold,
                            &query,
                            &locations,
                            &tables,
                            NCBI_OVERFLOW_CELL_CAPACITY,
                            None,
                            &[],
                            None,
                            passes,
                            stripe_fill,
                            Some((&spare, build == 1)),
                        )));
                        x_assert_same_flat_tables(&built, &serial);
                        if room && total > 0 {
                            assert_eq!(
                                built.hits.as_ptr() as usize,
                                given,
                                "round {round}, case {case}, build {build}: the spare's buffer"
                            );
                            reused_seen += 1;
                        }
                        let kept = spare.0.lock().unwrap();
                        assert!(
                            kept.0.is_empty() && kept.1.is_empty(),
                            "the build took them"
                        );
                        drop(kept);
                        spare.give(built);
                    }
                }
            }
        }
    }
    assert!(reused_seen > 0);
}

// LOSAT_X_BPLKPFLAT: the extra tasks run in the one-latch appends pass, each exactly once, at
// every thread count and with any number of chunks (none included), and not in the serial
// builder that a one-thread pool runs (`x_build_flat_in_pool`); the table does not change.
#[test]
fn extra_tasks_run_once_in_the_parallel_build_only() {
    let mut rng = XorShift(0xbb67_ae85_84ca_a73b);
    let lengths = [120, 7, 300];
    let (query, locations, first) = queries(&mut rng, &lengths, &FEW_LETTERS, 1);
    let tables = [first[0]..first[3]];
    let serial = flat(x_try_build_flat_compressed_lookup(
        5, 1100, &query, &locations, &tables,
    ));
    for threads in THREADS {
        for chunk_sizes in CHUNK_SIZES {
            for (count, stripe_fill) in [0, 1, 3, 64]
                .into_iter()
                .flat_map(|count| [(count, false), (count, true)])
            {
                let extra = Counted::new(count);
                let built = flat(Some(parallel(
                    threads,
                    1100,
                    &query,
                    &locations,
                    &tables,
                    NCBI_OVERFLOW_CELL_CAPACITY,
                    None,
                    chunk_sizes,
                    None,
                    XPasses::OneLatch(Some(&extra)),
                    stripe_fill,
                    None,
                )));
                x_assert_same_flat_tables(&built, &serial);
                extra.assert_each_ran(1, &format!("threads {threads}, {count} tasks"));
            }
        }
    }
    // No offsets: no chunk, the tasks still run once.
    for stripe_fill in [false, true] {
        let extra = Counted::new(5);
        let none = flat(Some(parallel(
            3,
            1100,
            &query,
            &[],
            &[],
            NCBI_OVERFLOW_CELL_CAPACITY,
            None,
            &[],
            None,
            XPasses::OneLatch(Some(&extra)),
            stripe_fill,
            None,
        )));
        assert!(none.hits.is_empty());
        extra.assert_each_ran(1, "no chunks");
    }
    // The entry of the search: the pool's build at 3 threads, the serial builder at one thread
    // (no task runs; the caller runs them), and nothing with the switch off.
    for (threads, mode, runs) in [
        (1, XLkpFlat::On, 0),
        (3, XLkpFlat::On, 1),
        (3, XLkpFlat::Off, 0),
    ] {
        let extra = Counted::new(9);
        let built = crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            Ok(super::super::x_build_flat_in_pool(
                5,
                1100,
                &query,
                &locations,
                &tables,
                pool,
                XFlatArgs {
                    mode,
                    spare: None,
                    extra: Some(&extra),
                    segments: None,
                    x_ahead: false,
                    appends: 0,
                },
            ))
        })
        .expect("pool");
        x_assert_same_flat_tables(&flat(built), &serial);
        extra.assert_each_ran(runs, &format!("threads {threads}, {mode:?}"));
    }
}

/// LOSAT_X_BPLKPAHEAD: the one-latch build ahead of its block (`x_ahead`) on `spare`.
#[allow(clippy::too_many_arguments)]
fn ahead(
    threads: usize,
    threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    spare: Option<&XFlatSpare>,
    poisoned: bool,
    stripe_fill: bool,
) -> Option<XFlatBuild> {
    if let (Some(spare), true) = (spare, poisoned) {
        poison(spare);
    }
    in_pool(threads, || {
        x_try_build_flat_parallel_with(
            5,
            threshold,
            query,
            locations,
            tables,
            XParallelCuts {
                overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
                max_appends: None,
                threads,
                chunk_sizes: &[],
                stripe_shift: None,
                wave_cells: None,
                appends_segment: X_APPENDS_SEGMENT,
                piece_appends: X_PIECE_APPENDS,
                segments: None,
                passes: XPasses::OneLatch(None),
                stripe_fill,
                spare,
                lent: None,
                x_ahead: true,
            },
        )
    })
}

/// The addresses of the spare's three vectors (`starts`, `hits`, `pv`).
fn spare_addresses(spare: &XFlatSpare) -> [usize; 3] {
    let kept = spare.0.lock().unwrap();
    let pv = spare.1.lock().unwrap();
    [
        kept.0.as_ptr() as usize,
        kept.1.as_ptr() as usize,
        pv.as_ptr() as usize,
    ]
}

/// Item SC2d: the parallel build's two passes run apart: the appends pass on a pool of
/// `threads.0` threads (`spare_a`: the spare it drops on `TooLarge`), its values moved to a pool of
/// `threads.1` threads for the stripes pass on `spare_b` (`x_ahead`: all of its vectors from
/// there).
#[allow(clippy::too_many_arguments)]
fn apart(
    threads: (usize, usize),
    threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    capacity: usize,
    max_appends: Option<usize>,
    passes: XPasses<'_>,
    stripe_fill: bool,
    spare_a: Option<&XFlatSpare>,
    spare_b: Option<&XFlatSpare>,
    x_ahead: bool,
) -> Option<XFlatBuild> {
    let stage = in_pool(threads.0, || {
        x_flat_appends_on(
            Blosum62Compressed,
            5,
            threshold,
            query,
            locations,
            tables,
            XParallelCuts {
                overflow_cell_capacity: capacity,
                max_appends,
                threads: threads.0,
                chunk_sizes: &[],
                stripe_shift: None,
                wave_cells: None,
                appends_segment: X_APPENDS_SEGMENT,
                piece_appends: X_PIECE_APPENDS,
                segments: None,
                passes,
                spare: spare_a,
                lent: None,
                x_ahead,
                stripe_fill,
            },
            |capacity| {
                x_try_build_flat_with_capacity(
                    5,
                    threshold,
                    query,
                    locations,
                    tables,
                    capacity,
                    X_COUNT_CHUNK,
                    0,
                )
            },
        )
    });
    match stage {
        XFlatStage::Built(built) => built,
        XFlatStage::Appends(appends) => in_pool(threads.1, || {
            x_flat_stripes_on(
                Blosum62Compressed,
                appends,
                query,
                locations,
                spare_b,
                None,
                x_ahead,
                None,
            )
        }),
    }
}

// Item SC2d: the appends pass and the stripes pass run apart (on pools of other sizes, the
// appends' values moved between them) give what the build gives at once: the table slot for slot,
// the refusal (the capacity test makes the words again) and `TooLarge`, with each pass form, the
// extra tasks run once in the appends pass, on fresh vectors, on a spare of garbage and (a build
// ahead of its block) on all of a spare's vectors. A spare of a build ahead without the room gives
// no table and keeps its vectors; a `TooLarge` build drops the appends pass's spare.
#[test]
fn the_two_passes_apart_build_what_the_build_at_once_builds() {
    let mut rng = XorShift(0x9b05_688c_2b3e_6c1f);
    let (mut refused_seen, mut built_seen, mut too_large_seen) = (0usize, 0usize, 0usize);
    for round in 0..6 {
        let threshold = [1100, 1600, 2000][round % 3];
        let letters: &[u8] = if round % 2 == 0 {
            &FEW_LETTERS
        } else {
            &LETTERS
        };
        let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 5 + rng.below(350)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 2);
        // Two tables, or one after segments of no table (its capacity test then reads the words).
        let half = first[lengths.len() / 2];
        let tables = if round % 3 == 0 {
            vec![first[0]..half, half..first[lengths.len()]]
        } else {
            vec![half..first[lengths.len()]]
        };
        let serial = flat(x_try_build_flat_compressed_lookup(
            5, threshold, &query, &locations, &tables,
        ));
        let total = serial.hits.len();
        let backbone_size = serial.starts.len() - 1;
        // The capacities: NCBI's, one cell (refusals where cells overflow); the bound: none, and
        // one below the appends.
        for (capacity, max_appends) in [
            (NCBI_OVERFLOW_CELL_CAPACITY, None),
            (1, None),
            (NCBI_OVERFLOW_CELL_CAPACITY, Some(total.saturating_sub(1))),
        ] {
            for (threads, (one_latch, stripe_fill)) in [
                ((2, 7), (false, false)),
                ((3, 2), (true, false)),
                ((7, 3), (true, true)),
            ] {
                let at_once = parallel(
                    threads.0,
                    threshold,
                    &query,
                    &locations,
                    &tables,
                    capacity,
                    max_appends,
                    &[],
                    None,
                    if one_latch {
                        XPasses::OneLatch(None)
                    } else {
                        XPasses::Rayon
                    },
                    stripe_fill,
                    None,
                );
                let extra = Counted::new(1 + rng.below(30));
                let passes = if one_latch {
                    XPasses::OneLatch(Some(&extra))
                } else {
                    XPasses::Rayon
                };
                // On a spare of garbage, given to the stripes pass.
                let spare = garbage_spare(&mut rng);
                give_room(&spare);
                poison(&spare);
                let built = apart(
                    threads,
                    threshold,
                    &query,
                    &locations,
                    &tables,
                    capacity,
                    max_appends,
                    passes,
                    stripe_fill,
                    Some(&spare),
                    Some(&spare),
                    false,
                )
                .expect("word size 5");
                let what = format!("round {round}, capacity {capacity}, bound {max_appends:?}, threads {threads:?}");
                match (&at_once, &built) {
                    (XFlatBuild::Table(want), XFlatBuild::Table(got)) => {
                        x_assert_same_flat_tables_as(&what, got, want);
                        built_seen += 1;
                    }
                    (XFlatBuild::Refused(want), XFlatBuild::Refused(got)) if want == got => {
                        refused_seen += 1
                    }
                    (XFlatBuild::TooLarge, XFlatBuild::TooLarge) => {
                        assert_eq!(spare_capacities(&spare), (0, 0), "{what}: dropped");
                        too_large_seen += 1;
                    }
                    _ => panic!("{what}: the outcomes differ"),
                }
                extra.assert_each_ran(usize::from(one_latch), &what);
                // A build ahead: all of the stripes pass's vectors from its spare (with garbage),
                // none dropped by the appends pass.
                if one_latch && capacity == NCBI_OVERFLOW_CELL_CAPACITY && max_appends.is_none() {
                    let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
                    let ahead_spare = XFlatSpare::default();
                    *ahead_spare.0.lock().unwrap() =
                        (vec![u32::MAX - 3; backbone_size + 1], vec![i32::MIN; total]);
                    *ahead_spare.1.lock().unwrap() = vec![u32::MAX; pv_words];
                    let given = spare_addresses(&ahead_spare);
                    let built = flat(apart(
                        threads,
                        threshold,
                        &query,
                        &locations,
                        &tables,
                        capacity,
                        None,
                        XPasses::OneLatch(None),
                        stripe_fill,
                        None,
                        Some(&ahead_spare),
                        true,
                    ));
                    x_assert_same_flat_tables_as(&what, &built, &serial);
                    let held = [
                        built.starts.as_ptr() as usize,
                        built.hits.as_ptr() as usize,
                        built.pv.as_ptr() as usize,
                    ];
                    assert_eq!(held, given, "{what}: the spare's vectors");
                    // One hit short of the room: no table, the spare keeps its vectors.
                    if total > 0 {
                        let short = XFlatSpare::default();
                        *short.0.lock().unwrap() = (
                            Vec::with_capacity(backbone_size + 1),
                            Vec::with_capacity(total - 1),
                        );
                        *short.1.lock().unwrap() = Vec::with_capacity(pv_words);
                        let kept = spare_addresses(&short);
                        let none = apart(
                            threads,
                            threshold,
                            &query,
                            &locations,
                            &tables,
                            capacity,
                            None,
                            XPasses::OneLatch(None),
                            stripe_fill,
                            None,
                            Some(&short),
                            true,
                        );
                        assert!(none.is_none(), "{what}: no table without the room");
                        assert_eq!(spare_addresses(&short), kept, "{what}: kept");
                    }
                }
            }
        }
    }
    assert!(refused_seen > 0 && built_seen > 0 && too_large_seen > 0);
}

// LOSAT_X_BPLKPAHEAD: a build ahead of its block takes `starts`, `hits` and the PV array from the
// spare that `x_make_room` made ready (garbage in them, the hits shorter than the table's) and
// allocates none of them: the table holds the spare's three buffers and equals the fresh serial
// build slot for slot, at 1, 2, 3 and 7 threads, poisoned or not. A spare one entry short of the
// room for any of the three (or no spare) gives no table and keeps its vectors untouched.
#[test]
fn builds_ahead_take_every_vector_from_the_spare() {
    let mut rng = XorShift(0x3c6e_f372_fe94_f82b);
    let mut built_seen = 0usize;
    for round in 0..4 {
        let (letters, threshold): (&[u8], i32) = match round % 2 {
            0 => (&FEW_LETTERS[..], 1100),
            _ => (&LETTERS[..], 1100),
        };
        let lengths: Vec<usize> = (0..2 + rng.below(4)).map(|_| 1 + rng.below(400)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 3);
        let tables = [first[0]..first[lengths.len()]];
        let serial = flat(x_try_build_flat_compressed_lookup(
            5, threshold, &query, &locations, &tables,
        ));
        let backbone_size = serial.starts.len() - 1;
        let total = serial.hits.len();
        let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
        for threads in THREADS {
            for (poison, stripe_fill) in
                [(false, false), (false, true), (true, false), (true, true)]
            {
                // The room, with garbage in it.
                let spare = XFlatSpare::default();
                let mut hits = Vec::with_capacity(total + 5);
                hits.extend((0..total / 2).map(|_| rng.below(1 << 20) as i32));
                *spare.0.lock().unwrap() = (vec![u32::MAX - 3; backbone_size + 1], hits);
                *spare.1.lock().unwrap() = vec![u32::MAX; pv_words];
                let block = XFlatSpare::default();
                assert_eq!(spare.x_make_room(&block, backbone_size, total, total), 0);
                let given = spare_addresses(&spare);
                let built = flat(ahead(
                    threads,
                    threshold,
                    &query,
                    &locations,
                    &tables,
                    Some(&spare),
                    poison,
                    stripe_fill,
                ));
                x_assert_same_flat_tables(&built, &serial);
                let held = [
                    built.starts.as_ptr() as usize,
                    built.hits.as_ptr() as usize,
                    built.pv.as_ptr() as usize,
                ];
                assert_eq!(held, given, "round {round}, threads {threads}");
                assert_eq!(spare_capacities(&spare), (0, 0));
                assert_eq!(spare.1.lock().unwrap().capacity(), 0);
                built_seen += 1;
                // One entry short of the room for each vector in turn, and no spare.
                for short in 0..3 {
                    if short == 1 && total == 0 {
                        continue;
                    }
                    let spare = XFlatSpare::default();
                    *spare.0.lock().unwrap() = (
                        Vec::with_capacity(backbone_size + usize::from(short != 0)),
                        Vec::with_capacity(total + usize::from(short != 1) - 1),
                    );
                    *spare.1.lock().unwrap() =
                        Vec::with_capacity(pv_words - usize::from(short == 2));
                    let kept = spare_addresses(&spare);
                    let none = ahead(
                        threads,
                        threshold,
                        &query,
                        &locations,
                        &tables,
                        Some(&spare),
                        poison,
                        stripe_fill,
                    );
                    assert!(none.is_none(), "round {round}, short {short}: no table");
                    assert_eq!(spare_addresses(&spare), kept, "the spare keeps them");
                }
                let none = ahead(
                    threads,
                    threshold,
                    &query,
                    &locations,
                    &tables,
                    None,
                    poison,
                    stripe_fill,
                );
                assert!(none.is_none());
            }
        }
    }
    assert!(built_seen > 0);
}

// LOSAT_X_BPLKPAHEAD: `x_make_room` keeps the spare's own vectors when they have the room, else
// moves the block spare's (the same buffers), replaces `starts` or the PV array short of its room
// by an empty one of the asked capacity and grows `hits` short of its room in place (item SC2d);
// it writes nothing into them.
#[test]
fn make_room_keeps_or_moves_the_vectors_and_allocates_only_what_is_short() {
    let backbone_size = 759_375;
    let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
    // Nothing kept anywhere: three new vectors.
    let ahead_spare = XFlatSpare::default();
    let block = XFlatSpare::default();
    assert_eq!(
        ahead_spare.x_make_room(&block, backbone_size, 1000, 1500),
        3
    );
    {
        let kept = ahead_spare.0.lock().unwrap();
        assert!(kept.0.is_empty() && kept.0.capacity() > backbone_size);
        assert!(kept.1.is_empty() && kept.1.capacity() >= 1500);
        assert!(ahead_spare.1.lock().unwrap().capacity() >= pv_words);
    }
    // Enough room: nothing allocated, the same buffers.
    let before = spare_addresses(&ahead_spare);
    assert_eq!(
        ahead_spare.x_make_room(&block, backbone_size, 1400, 3000),
        0
    );
    assert_eq!(spare_addresses(&ahead_spare), before);
    // `hits` short: only it is replaced (empty, with the grown capacity).
    assert_eq!(
        ahead_spare.x_make_room(&block, backbone_size, 1600, 3200),
        1
    );
    let after = spare_addresses(&ahead_spare);
    assert_eq!([after[0], after[2]], [before[0], before[2]]);
    assert!(ahead_spare.0.lock().unwrap().1.capacity() >= 3200);
    // An empty spare takes the block spare's vectors (the same buffers, contents kept).
    let empty = XFlatSpare::default();
    *block.0.lock().unwrap() = (vec![7; backbone_size + 1], vec![-7; 5000]);
    let moved = spare_addresses(&block);
    assert_eq!(empty.x_make_room(&block, backbone_size, 5000, 7500), 1);
    let taken = spare_addresses(&empty);
    assert_eq!([taken[0], taken[1]], [moved[0], moved[1]]);
    assert_eq!(spare_capacities(&block), (0, 0));
    assert_eq!(empty.0.lock().unwrap().1.len(), 5000);
}

// Item SC2d-2: the chunks' packed appends in the tasks' segments (`XSegmentOut`) of any size and
// in pieces of any size (`x_chunk_appends`): a task's pieces fill a segment and go on in a new one
// (a piece never straddles two), a piece of more appends than a segment holds gets one of its
// own, a chunk of more appends than a piece is cut into pieces, and the table is the serial one
// slot for slot, with rayon's passes and the one-latch passes (with and without the stripe fill),
// at 2, 3 and 7 threads, for BLOSUM62's words and another matrix's.
#[test]
fn the_appends_in_segments_and_pieces_of_any_size_build_the_serial_table() {
    let mut rng = XorShift(0x510e_527f_ade6_82d1);
    let (mut many_seen, mut whole_seen, mut cut_seen) = (0usize, 0usize, 0usize);
    for round in 0..4 {
        let matrix = [ScoringMatrix::Blosum62, ScoringMatrix::Blosum45][round % 2];
        let letters: &[u8] = if round < 2 { &FEW_LETTERS } else { &LETTERS };
        let lengths: Vec<usize> = (0..2 + rng.below(4)).map(|_| 20 + rng.below(300)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 2);
        let tables = vec![first[0]..first[lengths.len()]];
        let serial = flat(x_try_build_flat_of(
            matrix,
            5,
            1100,
            &query,
            &locations,
            &tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
            X_COUNT_CHUNK,
        ));
        for (segment, piece) in [
            (1, X_PIECE_APPENDS),
            (7, 1),
            (64, 5),
            (1000, 64),
            (X_APPENDS_SEGMENT, 300),
            (X_APPENDS_SEGMENT, X_PIECE_APPENDS),
        ] {
            for (threads, passes, stripe_fill) in [
                (2, XPasses::Rayon, false),
                (3, XPasses::OneLatch(None), false),
                (7, XPasses::OneLatch(None), true),
            ] {
                let what =
                    format!("round {round} {matrix:?} segment {segment} piece {piece} n{threads}");
                let cuts = XParallelCuts::<i32> {
                    overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
                    max_appends: None,
                    threads,
                    // Chunks of a few offsets, so that a task fills its segments with many.
                    chunk_sizes: &[3, 11, 29],
                    stripe_shift: None,
                    wave_cells: None,
                    appends_segment: segment,
                    piece_appends: piece,
                    segments: None,
                    passes,
                    spare: None,
                    lent: None,
                    x_ahead: false,
                    stripe_fill,
                };
                let stage = in_pool(threads, || {
                    x_flat_appends_on(matrix, 5, 1100, &query, &locations, &tables, cuts, |_| None)
                });
                let XFlatStage::Appends(appends) = stage else {
                    panic!("{what}: the appends pass ended the build")
                };
                // Each piece's range lies in its segment; a segment holds one piece larger
                // than `segment`, or pieces whose appends sum to at most `segment`. A piece
                // holds at most `piece` appends past its last offset's; the pieces are in
                // chunk order.
                let segments = &appends.appends.segments;
                let pieces = &appends.appends.pieces;
                let mut used = vec![0usize; segments.len()];
                for piece_appends in pieces {
                    assert!(piece_appends.len > 0, "{what}: an empty piece");
                    assert!(
                        piece_appends.at + piece_appends.len
                            <= segments[piece_appends.segment].len(),
                        "{what}: a piece past its segment"
                    );
                    used[piece_appends.segment] += piece_appends.len;
                    if piece_appends.len > segment {
                        assert_eq!(
                            segments[piece_appends.segment].len(),
                            piece_appends.len,
                            "{what}: a piece larger than a segment, in one of its own"
                        );
                        whole_seen += 1;
                    }
                }
                assert!(
                    pieces.windows(2).all(|w| w[0].chunk <= w[1].chunk),
                    "{what}: the pieces out of chunk order"
                );
                if pieces.windows(2).any(|w| w[0].chunk == w[1].chunk) {
                    cut_seen += 1;
                }
                for (index, (segment_entries, used)) in
                    segments.iter().map(Vec::len).zip(&used).enumerate()
                {
                    assert!(
                        *used <= segment_entries,
                        "{what}: segment {index} overfilled"
                    );
                }
                if segments.len() > threads {
                    many_seen += 1;
                }
                let built = in_pool(threads, || {
                    x_flat_stripes_on(matrix, appends, &query, &locations, None, None, false, None)
                });
                x_assert_same_flat_tables_as(&what, &flat(built), &serial);
            }
        }
    }
    assert!(
        many_seen > 0 && whole_seen > 0 && cut_seen > 0,
        "{many_seen} {whole_seen} {cut_seen}"
    );
}

// Item SC2d-2: guided chunk sizes (`XChunkSizes::Guided`): the chunks hold every offset of the
// count ranges once, in order, each within its range and `X_OFFSET_SPAN` of its first offset; a
// chunk takes the offsets left over `X_GUIDED_SHARE * threads` (at least `X_MIN_CHUNK_OFFSETS`)
// unless a range ends first, so the first chunks are the largest.
#[test]
fn guided_chunks_cover_the_offsets_in_order_and_shrink_to_the_floor() {
    let mut rng = XorShift(0x9b05_688c_2b3e_6c1e);
    for round in 0..6 {
        let lengths: Vec<usize> = (0..1 + rng.below(40))
            .map(|_| 1 + rng.below(3000))
            .collect();
        let (_, locations, first) = queries(&mut rng, &lengths, &LETTERS, round % 2);
        let half = first[lengths.len() / 2];
        let tables = if round % 2 == 0 {
            vec![first[0]..half, half..first[lengths.len()]]
        } else {
            vec![first[0]..first[lengths.len()]]
        };
        let ranges = x_count_ranges(&tables, locations.len());
        let word = 5i64;
        let expected: Vec<(usize, i64)> = ranges
            .iter()
            .enumerate()
            .flat_map(|(index, (range, _))| {
                locations[range.clone()]
                    .iter()
                    .flat_map(move |&(left, right)| {
                        (i64::from(left)..=i64::from(right) - word + 1).map(move |o| (index, o))
                    })
            })
            .collect();
        for threads in [1, 3, 8] {
            let sizes = XChunkSizes::Guided {
                threads,
                offsets: expected.len(),
            };
            let chunks = x_chunks(5, &ranges, &locations, sizes);
            let got: Vec<(usize, i64)> = chunks
                .iter()
                .flat_map(|chunk| {
                    chunk.pieces.iter().flat_map(move |&(a, b)| {
                        (i64::from(a)..=i64::from(b) - word + 1).map(move |o| (chunk.range, o))
                    })
                })
                .collect();
            assert_eq!(got, expected, "round {round} n{threads}: the offsets");
            let mut placed = 0usize;
            for (index, chunk) in chunks.iter().enumerate() {
                let count: usize = chunk
                    .pieces
                    .iter()
                    .map(|&(a, b)| (i64::from(b) - word + 2 - i64::from(a)) as usize)
                    .sum();
                let span = chunk
                    .pieces
                    .iter()
                    .map(|&(_, b)| i64::from(b) - word + 1 - i64::from(chunk.base))
                    .max()
                    .unwrap_or(0);
                assert!(span < X_OFFSET_SPAN, "round {round} n{threads}: span");
                let size = sizes.size(index, placed);
                assert!(count <= size, "round {round} n{threads} chunk {index}");
                // A chunk short of its size ends at its range's end or where the next offset
                // is past its offset span.
                let closed_early = chunks.get(index + 1).is_none_or(|next| {
                    next.range != chunk.range
                        || i64::from(next.base) >= i64::from(chunk.base) + X_OFFSET_SPAN
                });
                assert!(
                    count == size || closed_early,
                    "round {round} n{threads} chunk {index}: {count} of {size}"
                );
                placed += count;
            }
        }
    }
}

// Item SC2d-2: the segments a build keeps (`XSegmentPool`, at most as many as threads) may be
// the next build's (here filled with garbage first), and every build gives the serial table slot
// for slot.
#[test]
fn kept_segments_of_any_values_leave_the_table_the_serial_one() {
    let mut rng = XorShift(0x1f83_d9ab_fb41_bd6b);
    let pool = XSegmentPool::default();
    for round in 0..4 {
        let lengths: Vec<usize> = (0..2 + rng.below(4)).map(|_| 20 + rng.below(300)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, &LETTERS, round % 2);
        let tables = vec![first[0]..first[lengths.len()]];
        let serial = flat(x_try_build_flat_compressed_lookup(
            5, 1100, &query, &locations, &tables,
        ));
        let threads = 3;
        {
            let mut kept = pool.0.lock().unwrap();
            assert!(
                kept.len() <= threads,
                "round {round}: at most as many as threads"
            );
            for segment in kept.iter_mut() {
                segment.fill(u32::MAX);
            }
        }
        let stage = in_pool(threads, || {
            x_flat_appends_on(
                Blosum62Compressed,
                5,
                1100,
                &query,
                &locations,
                &tables,
                XParallelCuts::<i32> {
                    overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
                    max_appends: None,
                    threads,
                    chunk_sizes: &[3, 11, 29],
                    stripe_shift: None,
                    wave_cells: None,
                    appends_segment: 64,
                    piece_appends: X_PIECE_APPENDS,
                    segments: Some(&pool),
                    passes: XPasses::OneLatch(None),
                    spare: None,
                    lent: None,
                    x_ahead: false,
                    stripe_fill: true,
                },
                |_| None,
            )
        });
        let XFlatStage::Appends(appends) = stage else {
            panic!("round {round}: the appends pass ended the build")
        };
        // Which kept segments a build takes depends on the tasks' order (and a kept segment
        // shorter than a piece is not taken): `a_task_writes_into_kept_segments_first` checks
        // the order on one task. Here the kept ones hold u32::MAX and the table is the serial one.
        let built = in_pool(threads, || {
            x_flat_stripes_on(
                Blosum62Compressed,
                appends,
                &query,
                &locations,
                None,
                None,
                false,
                Some(&pool),
            )
        });
        x_assert_same_flat_tables_as(&format!("round {round}"), &flat(built), &serial);
        let kept = pool.0.lock().unwrap();
        assert!(
            (1..=threads).contains(&kept.len()),
            "round {round}: {} kept",
            kept.len()
        );
        assert!(
            kept.iter()
                .all(|segment| !segment.is_empty() && segment.capacity() == segment.len()),
            "round {round}: the kept segments are cut to their entries"
        );
    }
}

#[test]
fn a_task_writes_into_kept_segments_first() {
    let lengths = |pool: &XSegmentPool| -> Vec<usize> {
        pool.0.lock().unwrap().iter().map(Vec::len).collect()
    };
    let pool = XSegmentPool::default();
    // The empty segment is not kept.
    pool.keep(vec![vec![7u32; 10], Vec::new(), vec![7u32; 4]], 3);
    assert_eq!(lengths(&pool), [10, 4]);
    let kept: Vec<usize> = pool
        .0
        .lock()
        .unwrap()
        .iter()
        .map(|segment| segment.as_ptr() as usize)
        .collect();
    let registry = Mutex::new(Vec::new());
    {
        let mut out = XSegmentOut::new(&registry, 64, Some(&pool));
        let mut take = |n: usize, value: u32| {
            let (segment, at, range) = out.take(n);
            assert_eq!(range.len(), n);
            range.fill(value);
            (segment, at)
        };
        // The first kept segment of at least 6 entries (10), then 3 more in it.
        assert_eq!(take(6, 1), (0, 0));
        assert_eq!(take(3, 2), (0, 6));
        // 2 do not fit in the 1 left: the next kept one of at least 2 (4 entries).
        assert_eq!(take(2, 3), (1, 0));
        // 5 do not fit, and the pool is empty: a new segment of at least 64 entries.
        assert_eq!(take(5, 4), (2, 0));
    }
    let segments = registry.into_inner().unwrap();
    assert_eq!(
        segments.iter().map(Vec::len).collect::<Vec<_>>(),
        [9, 2, 5],
        "each segment cut to the entries taken"
    );
    assert_eq!(segments[0], [1, 1, 1, 1, 1, 1, 2, 2, 2]);
    assert_eq!(segments[1], [3, 3]);
    assert_eq!(segments[2], [4; 5]);
    assert_eq!(segments[0].as_ptr() as usize, kept[0]);
    assert_eq!(segments[1].as_ptr() as usize, kept[1]);
    assert!(lengths(&pool).is_empty());
    // At most 2 kept, each cut to its entries.
    pool.keep(segments, 2);
    assert_eq!(lengths(&pool), [9, 2]);
    assert!(pool
        .0
        .lock()
        .unwrap()
        .iter()
        .all(|segment| segment.capacity() == segment.len()));
    // A later build's segments replace the ones it did not take.
    pool.keep(vec![vec![5u32; 3]], 2);
    assert_eq!(lengths(&pool), [3]);

    // Item SC9i: between the waves of a build the segments are kept whole (at most 2): each at
    // the length it was made with, its blocks and entries the same.
    let mut made = vec![7u32; 12];
    made[..4].fill(8);
    made.truncate(4);
    let whole = (made.as_ptr() as usize, made.capacity());
    pool.keep_whole(vec![Vec::new(), made, vec![9u32; 5], vec![6u32; 2]], 2);
    let kept = pool.0.lock().unwrap();
    assert_eq!(kept.len(), 2);
    assert_eq!((kept[0].as_ptr() as usize, kept[0].len()), whole);
    assert_eq!(kept[0][..4], [8; 4]);
    assert_eq!(kept[0][4..], [7; 8]);
    assert_eq!(kept[1], [9; 5]);
}
