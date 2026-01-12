//! Offset pair storage.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-115

/// Query/subject offset pair (NCBI `BlastOffsetPair` equivalent for the hot path).
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct OffsetPair {
    pub q_off: u32,
    pub s_off: u32,
}
