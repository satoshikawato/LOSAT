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
//! `real` holds the real-number digits of every object stream; the Seq-annot model and the
//! encoders of item OF3 (outfmt-20261019) build on it.

pub mod real;
mod x87;

#[cfg(test)]
mod real_tests;
