//! BLAST XML2 and JSON (outfmt 15 and 16, and 13 and 14 as files; `CBlastOutput2`): the records of
//! NCBI's XML2 model beyond outfmt 5's, the type writers over `report::serial`'s encoders and
//! the four framings of a document (one document per query, joined for 15 and 16, written to
//! files for 13 and 14).
//!
//! The program adapters build the same per-query records as for outfmt 5
//! (`report::blast_xml`), which `Xml2Writer` writes through the `StructuredXml` trait. Design
//! and probes of NCBI BLAST+ 2.17.0: DESIGN_OF5 (`outfmt-20261019/S`).

pub mod files;
pub mod framing;
pub mod records;
pub mod writer;

pub use files::{base_name, FsSideFiles};
pub use framing::{NoSideFiles, RunSideFiles, SideFiles, Xml2Carry, Xml2Document, Xml2Kind};
pub use records::{subject_names, Xml2Run, Xml2Target};
pub use writer::Xml2Writer;

#[cfg(test)]
mod tests;
