/// `tests/tblastx_seed_compact.rs` counts the bytes that a compact seed record saves per init HSP
/// as `SAVED_PER_RECORD`, the size of an `InitHSP` minus the 28 bytes of the record; this test
/// ties that constant to the engine's struct.
#[cfg(target_pointer_width = "64")]
#[test]
fn init_hsp_size_is_the_one_the_seed_compact_test_counts() {
    assert_eq!(std::mem::size_of::<super::InitHSP>() - 28, 44);
}
