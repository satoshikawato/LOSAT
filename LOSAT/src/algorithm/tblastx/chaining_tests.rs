// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:98-100
// No NCBI counterpart: the layout of LOSAT's record. With the `Int4` values of NCBI in `u32`
// fields (and the subject frame index in `u8`), one ungapped HSP of the linking list takes
// 88 bytes on 64-bit targets instead of 152. On 32-bit targets (wasm32) `usize` was 4 bytes
// already, so nothing is pinned there.
#[cfg(target_pointer_width = "64")]
#[test]
fn ungapped_hit_is_88_bytes_on_64_bit_targets() {
    assert_eq!(std::mem::size_of::<super::UngappedHit>(), 88);
}
