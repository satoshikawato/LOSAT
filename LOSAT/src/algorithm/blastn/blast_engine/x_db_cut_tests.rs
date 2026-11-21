use super::*;

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

// The kernel runs once on every sequence and the kept results are merged once, in ordinal
// order.
#[test]
fn the_cut_searches_every_sequence_once_and_merges_in_ordinal_order() {
    let database = open("nucl_mv");
    let sequences = u64::from(database.db.num_seqs());
    assert_eq!(sequences, 12);
    let searched = std::sync::atomic::AtomicU64::new(0);
    let mut merged = Vec::new();
    crate::utils::threading::with_search_pool(2, "blastn", |pool| {
        x_search_blocks_cut(
            &database,
            pool,
            true,
            2,
            &|| 0u32,
            &|calls: &mut u32, ordinal: u32, _subject: BlastnSubjectRef<'_>| {
                *calls += 1;
                searched.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                ordinal.is_multiple_of(3).then_some(ordinal)
            },
            &mut |ordinal, kept| {
                assert_eq!(ordinal, kept);
                merged.push(ordinal);
            },
        )
    })
    .unwrap();
    assert_eq!(
        searched.load(std::sync::atomic::Ordering::Relaxed),
        sequences
    );
    assert_eq!(merged, [0, 3, 6, 9]);
}

// The cut: one thread keeps `BLOCK_BYTES`; more threads a quarter of the bytes per thread,
// between 4 KiB and 1 MiB.
#[test]
fn the_cut_is_a_quarter_of_the_bytes_per_thread_between_4_kib_and_1_mib() {
    let x_cut_block_bytes = |bytes, threads| x_cut_block_bytes(bytes, threads, BLOCK_BYTES);
    for bytes in [0u64, 1, 812, 10_300_000, 64 << 20, u64::MAX] {
        assert_eq!(x_cut_block_bytes(bytes, 0), BLOCK_BYTES, "{bytes}");
        assert_eq!(x_cut_block_bytes(bytes, 1), BLOCK_BYTES, "{bytes}");
    }
    // 28,036 16S rRNA genes (10.3 MB) at 16 threads: 64 blocks of 160,937 bytes; their first
    // 2,000 (0.72 MB) at 16 threads: blocks of 11,250 bytes.
    assert_eq!(x_cut_block_bytes(10_300_000, 16), 160_937);
    assert_eq!(x_cut_block_bytes(720_000, 16), 11_250);
    assert_eq!(x_cut_block_bytes(720_000, 2), 90_000);
    // Small databases: the floor.
    assert_eq!(x_cut_block_bytes(0, 2), X_MIN_BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(8 * 4096 - 1, 2), X_MIN_BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(8 * 4096 + 8, 2), 4097);
    // From 4 x threads MiB on: `BLOCK_BYTES`.
    assert_eq!(x_cut_block_bytes((64 << 20) - 64, 16), BLOCK_BYTES - 1);
    assert_eq!(x_cut_block_bytes(64 << 20, 16), BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(u64::MAX, 7), BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(u64::MAX, usize::MAX), X_MIN_BLOCK_BYTES);
    // Never above one thread's blocks, never below the floor; smaller with more threads.
    for bytes in [0u64, 5_000, 300_000, 10_300_000, 40 << 20, 200 << 20] {
        let mut previous = BLOCK_BYTES;
        for threads in 1..=130 {
            let cut = x_cut_block_bytes(bytes, threads);
            assert!((X_MIN_BLOCK_BYTES..=BLOCK_BYTES).contains(&cut));
            assert!(cut <= previous, "{bytes} bytes, {threads} threads");
            previous = cut;
        }
    }
}

// LOSAT_X_BNDBSMALL: the largest block only. One thread keeps `BLOCK_BYTES`, a smaller cut stays.
#[test]
fn small_blocks_only_lower_the_largest_cut() {
    let small = X_BNDBSMALL_BLOCK_BYTES;
    assert_eq!(x_cut_block_bytes(u64::MAX, 1, small), BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(344_000_000, 16, small), small);
    assert_eq!(x_cut_block_bytes(10_300_000, 16, small), 160_937);
    assert_eq!(x_cut_block_bytes(0, 2, small), X_MIN_BLOCK_BYTES);
}
