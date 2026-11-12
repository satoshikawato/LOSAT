//! NCBI's `objtools/seqmasks_io` (`c++/src/objtools/seqmasks_io/`): the input and output of
//! the masking programs (dustmasker, segmasker, windowmasker).
#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "read by the masking programs, which land after this module (line RG, K1)"
    )
)]

pub(crate) mod fasta_reader;
