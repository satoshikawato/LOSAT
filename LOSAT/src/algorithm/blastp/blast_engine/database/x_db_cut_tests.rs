use super::*;

// The cut: one thread keeps `BLOCK_BYTES`; more threads a quarter of the bytes per thread,
// between 4 KiB and 1 MiB.
#[test]
fn the_cut_is_a_quarter_of_the_bytes_per_thread_between_4_kib_and_1_mib() {
    for bytes in [0u64, 1, 812, 2_920_000, 64 << 20, u64::MAX] {
        assert_eq!(x_cut_block_bytes(bytes, 0), BLOCK_BYTES, "{bytes}");
        assert_eq!(x_cut_block_bytes(bytes, 1), BLOCK_BYTES, "{bytes}");
    }
    // BlastDBSpecies6 (2.92 MB) at 16 threads: 64 blocks of 45,625 bytes.
    assert_eq!(x_cut_block_bytes(2_920_000, 16), 45_625);
    assert_eq!(x_cut_block_bytes(2_920_000, 2), 365_000);
    // Small databases: the floor.
    assert_eq!(x_cut_block_bytes(0, 2), X_MIN_BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(812, 7), X_MIN_BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(8 * 4096 - 1, 2), X_MIN_BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(8 * 4096 + 8, 2), 4097);
    // From 4 x threads MiB on: `BLOCK_BYTES`.
    assert_eq!(x_cut_block_bytes((64 << 20) - 64, 16), BLOCK_BYTES - 1);
    assert_eq!(x_cut_block_bytes(64 << 20, 16), BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(u64::MAX, 7), BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(u64::MAX, usize::MAX), X_MIN_BLOCK_BYTES);
    assert_eq!(x_cut_block_bytes(1 << 40, usize::MAX), X_MIN_BLOCK_BYTES);
    // Never above one thread's blocks, never below the floor; smaller with more threads.
    for bytes in [0u64, 5_000, 300_000, 2_920_000, 40 << 20, 200 << 20] {
        let mut previous = BLOCK_BYTES;
        for threads in 1..=130 {
            let cut = x_cut_block_bytes(bytes, threads);
            assert!((X_MIN_BLOCK_BYTES..=BLOCK_BYTES).contains(&cut));
            assert!(cut <= previous, "{bytes} bytes, {threads} threads");
            previous = cut;
        }
    }
}
