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

pub use records::{HitTaxonomy, QueryMask, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats};
pub use subjects::{hit_ranges, SubjectNameCache, XmlSubjects, XmlTitles};
pub use writer::{written_queries, BlastXmlWriter, DocumentPosts, StructuredXml};

/// An exception of NCBI's formatters raised while one query's structured report is made, which
/// NCBI's report function catches itself (BLAST XML's `s_BlastXMLAddIteration`, the stream forms
/// of `BlastXML2_FormatReport` and `BlastJSON_FormatReport`): the document writes what NCBI wrote
/// before it (`StructuredXml::query`'s output) and the program posts its message
/// (`StructuredXml::caught_message`) and goes on. Raised by BLASTX's HSPs of length 0
/// (`algorithm::blastx::empty_rows`) and by the `Hit` of a database sequence whose lookup throws
/// where the `Hit` loads it (`subjects::HitNames::names` of `HitNames::Throws`, as
/// `report::db_subject::HitLoads` decides it; XML2 and JSON only).
#[derive(Debug)]
pub struct CaughtException(pub crate::ncbi_exception::NcbiException);

impl std::fmt::Display for CaughtException {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(self.0.get_msg()))
    }
}

impl std::error::Error for CaughtException {}

/// `BlastOutput_reference` of a Gapped BLAST search (every task but megablast): NCBI's text
/// (`CReference::GetString` of `eGappedBlast`, `report::reference`), whose `&auml;` is escaped
/// again when written.
pub const REFERENCE_GAPPED_BLAST: &str = super::reference::Publication::GappedBlast.string();

/// `BlastOutput_reference` of megablast (`CReference::GetString` of `eMegaBlast`).
pub const REFERENCE_MEGABLAST: &str = super::reference::Publication::MegaBlast.string();

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:652-671
/// `BlastOutput_reference` of psiblast, whose task is `ePSIBlast` in every round
/// (`CReference::GetString` of `eCompBasedStats`).
pub const REFERENCE_COMP_BASED_STATS: &str = super::reference::Publication::CompBasedStats.string();

#[cfg(test)]
mod tests;
