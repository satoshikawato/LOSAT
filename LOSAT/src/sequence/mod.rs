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
/// (`blastdbcmd::fasta`).
pub(crate) fn write_fasta_lines(letters: &[u8], width: u32, out: &mut Vec<u8>) {
    let mut lines = FastaLines::new(width);
    lines.write(letters, out);
    lines.end(out);
}

/// The lines of [`write_fasta_lines`] written as the letters come, a part at a time: the
/// queries and subjects of a search strategy (`report::archive::read::queries`), which go on to
/// a file a block at a time.
pub(crate) struct FastaLines {
    width: usize,
    /// The letters of the line being written.
    column: usize,
}

impl FastaLines {
    pub(crate) fn new(width: u32) -> Self {
        Self {
            width: width.max(1) as usize,
            column: 0,
        }
    }

    /// The next `letters` of the sequence, a line feed after each full line.
    pub(crate) fn write(&mut self, mut letters: &[u8], out: &mut Vec<u8>) {
        while !letters.is_empty() {
            let taken = (self.width - self.column).min(letters.len());
            out.extend_from_slice(&letters[..taken]);
            letters = &letters[taken..];
            self.column += taken;
            if self.column == self.width {
                out.push(b'\n');
                self.column = 0;
            }
        }
    }

    /// The end of the sequence: the line feed of a last line that is not full.
    pub(crate) fn end(self, out: &mut Vec<u8>) {
        if self.column > 0 {
            out.push(b'\n');
        }
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
