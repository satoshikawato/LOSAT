//! Sequence representation and manipulation utilities
//!
//! This module provides efficient sequence representations for BLAST operations,
//! including the 2-bit packed nucleotide format used by NCBI BLAST.

pub mod packed_letters;
pub mod packed_subject;
pub(crate) mod stdaa_letters;
pub(crate) mod title;

// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3295-3317
/// The sequence lines of `CFastaOstream::x_WriteSequence`: lines of `width` letters, the
/// last one shorter, each ended by a line feed (none for no letters): blastdbcmd's FASTA
/// (`blastdbcmd::fasta`) and the queries of a search strategy (`report::archive::read::queries`).
pub(crate) fn write_fasta_lines(letters: &[u8], width: u32, out: &mut Vec<u8>) {
    for line in letters.chunks(width.max(1) as usize) {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
