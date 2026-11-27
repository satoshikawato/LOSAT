use super::super::{search_rounds, BlastnDatabase};
use super::*;
use std::sync::atomic::AtomicU64;

fn open(name: &str) -> BlastnDatabase {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/db")
        .join(name);
    BlastnDatabase::open(
        path.to_str().expect("fixture path is UTF-8"),
        &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
    )
    .unwrap_or_else(|error| panic!("{name}: {error:#}"))
}

// LOSAT_X_BNDBSMALL: the window keeps its bytes.
#[test]
fn the_wide_window_keeps_its_bytes() {
    assert_eq!(x_wide_blocks_per_thread(BLOCK_BYTES), 16);
    assert_eq!(x_wide_blocks_per_thread(BLOCK_BYTES - 1), 16);
    assert_eq!(x_wide_blocks_per_thread(256 << 10), 64);
    assert_eq!(x_wide_blocks_per_thread(1), 16 << 20);
    assert_eq!(x_wide_blocks_per_thread(0), 16 << 20);
}

// Uneven work per sequence, so that the blocks end out of order.
// The bytes that wait are the results' lists and the list of the block's kept results.
#[test]
fn the_kept_bytes_are_the_capacities_of_the_lists() {
    let mut kept: Vec<(u32, Vec<u64>)> = Vec::with_capacity(3);
    kept.push((0, Vec::with_capacity(100)));
    kept.push((1, Vec::with_capacity(1000)));
    let outer = kept.capacity() * std::mem::size_of::<(u32, Vec<u64>)>();
    let inner = (kept[0].1.capacity() + kept[1].1.capacity()) * std::mem::size_of::<u64>();
    assert_eq!(kept_bytes(&kept), outer + inner);
}

fn work(ordinal: u32) -> u64 {
    let rounds = if ordinal % 5 == 1 { 200_000 } else { 100 };
    (0..rounds).fold(u64::from(ordinal), |acc, x| {
        acc.wrapping_mul(31).wrapping_add(x)
    })
}

// Every sequence is searched once, the kept results reach `merge` once each in ordinal
// order, and a thread makes one scratch at most: at 1 to 8 threads, blocks of one sequence
// (1 byte) to the whole database.
#[test]
fn the_window_merges_in_ordinal_order_once() {
    let database = open("nucl_mv");
    let db = &database.db;
    let sequences = db.num_seqs();
    assert_eq!(sequences, 12);
    let want: Vec<(u32, u64)> = (0..sequences)
        .filter(|ordinal| ordinal % 3 != 2)
        .map(|ordinal| (ordinal, work(ordinal)))
        .collect();
    for threads in 1..=8 {
        for block_bytes in [1u64, 40, 1 << 20] {
            let searched = AtomicU64::new(0);
            let scratches = AtomicU64::new(0);
            let mut merged: Vec<(u32, u64)> = Vec::new();
            crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                x_search_blocks_stream(
                    &database,
                    pool,
                    block_bytes,
                    &|| {
                        scratches.fetch_add(1, Ordering::Relaxed);
                        Vec::<u32>::new()
                    },
                    &|seen: &mut Vec<u32>, ordinal: u32, subject: BlastnSubjectRef<'_>| {
                        assert!(subject.length > 0);
                        // a scratch sees its sequences in increasing order
                        assert!(seen.last().is_none_or(|&last| last < ordinal));
                        seen.push(ordinal);
                        searched.fetch_add(1, Ordering::Relaxed);
                        (ordinal % 3 != 2).then(|| work(ordinal))
                    },
                    &mut |ordinal, kept| merged.push((ordinal, kept)),
                )
            })
            .unwrap();
            let case = format!("threads {threads} block_bytes {block_bytes}");
            assert_eq!(merged, want, "{case}");
            assert_eq!(
                searched.load(Ordering::Relaxed),
                u64::from(sequences),
                "{case}"
            );
            assert!(
                (1..=threads as u64).contains(&scratches.load(Ordering::Relaxed)),
                "{case}: {} scratches",
                scratches.load(Ordering::Relaxed)
            );
        }
    }
}

// The list of blocks is the order `next_block` gives: consecutive, every sequence once.
#[test]
fn the_blocks_cover_the_database_in_order() {
    let database = open("nucl_mv");
    let db = &database.db;
    for block_bytes in [1u64, 40, 100, 1 << 20] {
        let blocks = blocks_of(db, 0, block_bytes, usize::MAX);
        let mut next = 0u32;
        for range in &blocks {
            assert_eq!(range.first_ordinal, next, "{block_bytes}");
            assert!(range.end_ordinal() > next, "{block_bytes}");
            next = range.end_ordinal();
        }
        assert_eq!(next, db.num_seqs(), "{block_bytes}");
    }
    assert_eq!(blocks_of(db, 0, 1, usize::MAX).len(), 12);
    // Segments: up to `max` blocks from the given ordinal id on.
    let first = blocks_of(db, 0, 1, 5);
    assert_eq!(first.len(), 5);
    let rest = blocks_of(db, first[4].end_ordinal(), 1, usize::MAX);
    assert_eq!(rest.len(), 7);
    assert_eq!(rest[0].first_ordinal, 5);
}

// Segments of 1, 5 and 64 blocks: every kept result once, in ordinal order.
#[test]
fn segments_merge_in_ordinal_order_once() {
    let database = open("nucl_mv");
    let db = &database.db;
    let want: Vec<u32> = (0..db.num_seqs()).filter(|o| o % 2 == 0).collect();
    for threads in [2usize, 5] {
        for segment in [1usize, 5, 64] {
            let mut merged = Vec::new();
            let result = crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                Ok(stream_blocks(
                    &database,
                    pool,
                    1,
                    segment,
                    &|| (),
                    &|_: &mut (), ordinal: u32, _: BlastnSubjectRef<'_>| {
                        ordinal.is_multiple_of(2).then_some(ordinal)
                    },
                    &mut |ordinal, kept| {
                        assert_eq!(ordinal, kept);
                        merged.push(ordinal);
                    },
                ))
            })
            .unwrap();
            result.unwrap();
            assert_eq!(merged, want, "threads {threads} segment {segment}");
        }
    }
}

// A volume in the middle of the database cut short: the window returns the rounds' error (the
// first failing block in ordinal order) after merging the same blocks before it, and nothing
// after it, at 2 to 8 threads with blocks of one sequence, of 40 bytes and of 1 MiB.
#[test]
fn a_failing_block_in_the_middle_stops_the_window_as_it_stops_the_rounds() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-bndbstream-middle-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db");
    for entry in std::fs::read_dir(&fixtures).unwrap() {
        let entry = entry.unwrap();
        let name = entry.file_name().into_string().unwrap();
        if name.starts_with("nucl_mv.") {
            std::fs::copy(entry.path(), dir.join(&name)).unwrap();
        }
    }
    let nsq = dir.join("nucl_mv.02.nsq");
    let bytes = std::fs::read(&nsq).unwrap();
    std::fs::write(&nsq, &bytes[..bytes.len() / 2]).unwrap();
    let database = BlastnDatabase::open(
        dir.join("nucl_mv").to_str().unwrap(),
        &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
    )
    .unwrap();
    let db = &database.db;
    let kernel = |_: &mut (), ordinal: u32, _: BlastnSubjectRef<'_>| Some(ordinal);
    for threads in 2..=8 {
        for block_bytes in [1u64, 40, 1 << 20] {
            let (mut rounds, mut window) = (Vec::new(), Vec::new());
            let (rounds_result, window_result) =
                crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                    let rounds_result = search_rounds(
                        &database,
                        pool,
                        block_bytes,
                        2,
                        &|| (),
                        &kernel,
                        &mut |o, k| {
                            assert_eq!(o, k);
                            rounds.push(o);
                        },
                    );
                    let window_result = stream_blocks(
                        &database,
                        pool,
                        block_bytes,
                        3,
                        &|| (),
                        &kernel,
                        &mut |o, k| {
                            assert_eq!(o, k);
                            window.push(o);
                        },
                    );
                    Ok((rounds_result, window_result))
                })
                .unwrap();
            let case = format!("threads {threads} block_bytes {block_bytes}");
            let rounds_error = format!("{:#}", rounds_result.expect_err(&case));
            let window_error = format!("{:#}", window_result.expect_err(&case));
            assert_eq!(window_error, rounds_error, "{case}");
            assert_eq!(window, rounds, "{case}");
            assert!(window.len() < db.num_seqs() as usize, "{case}");
            assert_eq!(
                window,
                (0..window.len() as u32).collect::<Vec<_>>(),
                "{case}"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}
