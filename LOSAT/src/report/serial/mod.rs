//! NCBI serial-library output (the object streams behind outfmt 5, 8, 9, 11, 12 and 13-16).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:694-705
//! ```c++
//!             if(m_FormatType == CFormattingArgs::eAsnText)
//!                 m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//!             else
//!                 m_Outfile << MSerial_Json << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//!         ...
//!             m_Outfile << MSerial_AsnBinary <<
//!                 *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//! ```
//!
//! `real` holds the real-number digits of every object stream; `encoder` the object-stream
//! calls that the type writers make, with one implementation per format (`asn_text`,
//! `json`, `asn_ber`); `seqalign` the Seq-annot of Seq-aligns of outfmt 8, 9 and 12 (`write_query`).

pub mod asn_ber;
pub mod asn_text;
pub mod encoder;
pub mod json;
pub mod real;
pub mod seqalign;
mod x87;

pub use seqalign::{write_query, SeqAlignError, SerialFormat};

#[cfg(test)]
mod asn_ber_tests;
#[cfg(test)]
mod real_tests;
#[cfg(test)]
mod seqalign_tests;
