//! The input and output of NCBI's masking tools (`c++/src/objtools/seqmasks_io/`): the mask
//! writers of `dustmasker` (and later `segmasker`, `windowmasker`), one for each `-outfmt`.
//!
//! A tool hands each sequence and its masked intervals to one [`writer::MaskWriter`], then
//! finishes it (the `maskinfo_*` formats write their one object at the end).

pub(crate) mod writer;
