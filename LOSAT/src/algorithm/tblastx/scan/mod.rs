//! Scan support for TBLASTX (offset pair storage).
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-115

mod offset_pairs;

pub use offset_pairs::OffsetPair;
