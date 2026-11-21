//! Offset pair storage.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:99-115

/// Query/subject offset pair (NCBI `BlastOffsetPair` equivalent for the hot path).
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct OffsetPair {
    pub q_off: u32,
    pub s_off: u32,
}
