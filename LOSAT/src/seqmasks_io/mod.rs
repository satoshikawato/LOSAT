//! NCBI's `objtools/seqmasks_io` (`c++/src/objtools/seqmasks_io/`): the input and output of
//! the masking programs (dustmasker; segmasker and windowmasker later).
//!
//! `fasta_reader` reads FASTA input (`CMaskFastaReader`). A tool hands each sequence and its
//! masked intervals to one [`writer::MaskWriter`], then finishes it (the `maskinfo_*` formats
//! write their one object at the end).

pub(crate) mod fasta_reader;
pub(crate) mod writer;
