//! Sequence representation and manipulation utilities
//!
//! This module provides efficient sequence representations for BLAST operations,
//! including the 2-bit packed nucleotide format used by NCBI BLAST.

pub mod packed_letters;
pub mod packed_subject;
pub(crate) mod stdaa_letters;

// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3295-3317
// ```c++
//             while ( count >= rem_line ) {
//                 x_WriteBuffer(ptr, rem_line);
//                 if ( !(current_state & eHardMask) ) {
//                     ptr += rem_line;
//                 }
//                 count -= rem_line;
//                 m_Out << '\n';
//                 rem_line = m_Width;
//             }
//             if ( count > 0 ) {
//                 x_WriteBuffer(ptr, count);
//                 rem_line -= count;
//             }
//             it.SetPos(new_pos);
//         }
//     }
//     if ( rem_line < m_Width ) {
//         m_Out << '\n';
//     }
// ```
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
mod tests {
    use super::*;

    #[test]
    fn sequence_lines_are_the_width_with_a_shorter_last() {
        let mut out = Vec::new();
        write_fasta_lines(b"ACGTACGT", 3, &mut out);
        assert_eq!(out, b"ACG\nTAC\nGT\n");
        out.clear();
        write_fasta_lines(b"ACGTAC", 3, &mut out);
        assert_eq!(out, b"ACG\nTAC\n");
    }
}
