//! NCBI's `objtools/seqmasks_io` (`c++/src/objtools/seqmasks_io/`): the input and output of
//! the masking programs (dustmasker, segmasker; windowmasker later).
//!
//! `fasta_reader` reads FASTA input (`CMaskFastaReader`), `bdb_reader` a BLAST database
//! (`CMaskBDBReader`); `run` holds the parts of the programs' runs around them. A tool hands
//! each sequence and its masked intervals to one
//! [`MaskWriter`](crate::seqmasks_io::writer::MaskWriter), then finishes it (the `maskinfo_*`
//! formats write their one object at the end).

pub(crate) mod bdb_reader;
pub(crate) mod fasta_reader;
pub(crate) mod run;
pub(crate) mod writer;
