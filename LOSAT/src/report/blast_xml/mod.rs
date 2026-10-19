//! BLAST XML (outfmt 5, `CBlastOutput`): the records of NCBI's XML model, the value rules of its
//! writer and the incremental document, over `report::serial::xml::XmlEncoder`.
//!
//! The program adapters (`algorithm/<program>/blast_xml.rs`) build `XmlQuery`, `XmlHit` and
//! `XmlHsp` from their final hit lists, per query, into buffers reused from query to query.
//! Design and probes of NCBI BLAST+ 2.17.0: DESIGN_OF4 (`outfmt-20261019/S`).

pub mod filter;
pub mod hsp;
pub mod matrix;
pub mod names;
pub mod records;
pub mod subjects;
pub mod writer;

pub use records::{QueryMask, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats};
pub use subjects::XmlSubjects;
pub use writer::{BlastXmlWriter, HitSink, HspSink};

/// `BlastOutput_reference` of a Gapped BLAST search (every task but megablast); the
/// `&auml;` is NCBI's text (api/version.cpp:47-51) and is escaped again when written.
pub const REFERENCE_GAPPED_BLAST: &str = "Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of protein database search programs\", Nucleic Acids Res. 25:3389-3402.";

/// `BlastOutput_reference` of megablast (`CReference::eMegaBlast`).
pub const REFERENCE_MEGABLAST: &str = "Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), \"A greedy algorithm for aligning DNA sequences\", J Comput Biol 2000; 7(1-2):203-14.";

#[cfg(test)]
mod tests;
