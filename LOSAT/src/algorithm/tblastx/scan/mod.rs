//! Scan support for TBLASTX (offset pair storage).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:99-115

mod offset_pairs;

pub use offset_pairs::OffsetPair;
