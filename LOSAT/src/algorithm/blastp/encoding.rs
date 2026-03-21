//! Protein sequence encoding helpers for BLASTP

use crate::algorithm::tblastx::translation::QueryFrame;
use crate::utils::matrix::{aa_char_to_ncbistdaa, ncbistdaa};
use crate::utils::seg::{SegMasker, SegParams};

/// Encoded protein sequence in NCBISTDAA with leading and trailing sentinels.
#[derive(Debug, Clone)]
pub struct EncodedProtein {
    pub aa_seq: Vec<u8>,
    pub aa_len: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:115-121
// ```c
// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
// '-','A','B','C','D','E','F','G','H','I','K','L','M',
// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
// 'O', 'J'};
// const Uint1 kProtSentinel = NULLB;
// ```
pub fn ncbistdaa_to_ascii(residue: u8) -> char {
    const TABLE: [char; 28] = [
        '-', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
        'T', 'V', 'W', 'X', 'Y', 'Z', 'U', '*', 'O', 'J',
    ];
    TABLE.get(residue as usize).copied().unwrap_or('X')
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-127
// ```c
// (*seq_blk)->sequence_start = (Uint1 *) buffer;
// /* The first byte is a sentinel byte. */
// (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
// (*seq_blk)->sequence_start_nomask = (*seq_blk)->sequence_start;
// (*seq_blk)->sequence_nomask = (*seq_blk)->sequence;
// (*seq_blk)->length = length;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:110-121
// ```c
// ... AMINOACID_TO_NCBISTDAA ...
// const Uint1 kProtSentinel = NULLB;
// ```
pub fn encode_protein_sequence(seq: &[u8]) -> EncodedProtein {
    let mut aa_seq = Vec::with_capacity(seq.len() + 2);
    aa_seq.push(0);
    for &residue in seq {
        aa_seq.push(aa_char_to_ncbistdaa(residue.to_ascii_uppercase()));
    }
    aa_seq.push(0);
    EncodedProtein {
        aa_seq,
        aa_len: seq.len(),
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:382-385
// ```c
// if (Blast_QueryIsProtein(program_number))
//     BLAST_ScoreSetAmbigRes(sbp, 'X');
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:337-370
// ```c
// BlastSetUp_Filter(..., query_blk, ...);
// query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
// query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
// ```
pub fn encode_protein_query_frame_with_seg(seq: &[u8], seg: Option<&SegParams>) -> QueryFrame {
    let encoded = encode_protein_sequence(seq);
    let mut aa_seq = encoded.aa_seq.clone();
    let mut aa_seq_nomask = None;
    let mut seg_masks = Vec::new();

    if let Some(params) = seg {
        let masker = SegMasker::with_params(params);
        let intervals = masker.mask_sequence(&aa_seq[1..aa_seq.len() - 1]);
        if !intervals.is_empty() {
            aa_seq_nomask = Some(aa_seq.clone());
            for interval in &intervals {
                seg_masks.push((interval.start, interval.end));
                for pos in interval.start..interval.end {
                    aa_seq[pos + 1] = ncbistdaa::X;
                }
            }
        }
    }

    QueryFrame {
        frame: 0,
        aa_seq,
        aa_seq_nomask,
        aa_len: encoded.aa_len,
        orig_len: encoded.aa_len,
        seg_masks,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-315
// ```c
// query_offset = 0;
// for (index = 0; index < query_info->last_context; ++index) {
//     query_info->contexts[index].query_offset = query_offset;
//     query_offset += query_info->contexts[index].query_length + 1;
// }
// ```
pub fn encode_protein_query_frame(seq: &[u8]) -> QueryFrame {
    encode_protein_query_frame_with_seg(seq, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-127
    // ```c
    // (*seq_blk)->sequence_start = (Uint1 *) buffer;
    // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
    // (*seq_blk)->length = length;
    // ```
    #[test]
    fn test_encode_protein_sequence_adds_sentinels() {
        let encoded = encode_protein_sequence(b"ACD");
        assert_eq!(encoded.aa_len, 3);
        assert_eq!(encoded.aa_seq.len(), 5);
        assert_eq!(encoded.aa_seq[0], 0);
        assert_eq!(encoded.aa_seq[4], 0);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:115-121
    // ```c
    // const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
    // '-','A','B','C','D','E','F','G','H','I','K','L','M',
    // 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
    // 'O', 'J'};
    // ```
    #[test]
    fn test_ncbistdaa_to_ascii_roundtrip() {
        let encoded = encode_protein_sequence(b"Acd");
        let letters: String = encoded.aa_seq[1..encoded.aa_seq.len() - 1]
            .iter()
            .map(|&b| ncbistdaa_to_ascii(b))
            .collect();
        assert_eq!(letters, "ACD");
    }
}
