//! Protein sequence encoding helpers for BLASTP

use crate::algorithm::blastx::query_setup::combine_masks;
use crate::algorithm::tblastn::archive::lowercase_runs;
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
// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:894-899
// ```c
//     for (i = 0; i < sv.size(); i++) {
//         // Change unsupported residues to X
//         if (sv[i] == AMINOACID_TO_NCBISTDAA[(int)'O']) {
//             replaced_residues.push_back(i);
//             *buf_var++ = AMINOACID_TO_NCBISTDAA[(int)'X'];
// ```
// The search and its statistics read pyrrolysine (O) as X; the reports show the input.
pub fn encode_protein_sequence(seq: &[u8]) -> EncodedProtein {
    let mut aa_seq = Vec::with_capacity(seq.len() + 2);
    aa_seq.push(0);
    for &residue in seq {
        let residue = match residue.to_ascii_uppercase() {
            b'O' => b'X',
            upper => upper,
        };
        aa_seq.push(aa_char_to_ncbistdaa(residue));
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
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-315
// ```c
// query_offset = 0;
// for (index = 0; index < query_info->last_context; ++index) {
//     query_info->contexts[index].query_offset = query_offset;
//     query_offset += query_info->contexts[index].query_length + 1;
// }
// ```
pub fn encode_protein_query_frame_with_seg(seq: &[u8], seg: Option<&SegParams>) -> QueryFrame {
    encode_query_frame::<false>(seq, seg, false)
}

/// As `encode_protein_query_frame_with_seg`, and with LOSAT_TIMING the encoding and the SEG are
/// added to the common line's `encode` and `seg` (utils/xtiming.rs). For the drivers' query
/// preparation on the calling thread only, not for per-subject work.
pub fn encode_protein_query_frame_with_seg_timed(
    seq: &[u8],
    seg: Option<&SegParams>,
) -> QueryFrame {
    encode_query_frame::<false>(seq, seg, true)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
// ```c
//     if (!mask_at_hash) {
//         BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc,
//                              program_number);
//     }
// ```
/// The query frame of a BLASTP search with `-soft_masking true` (`mask_at_hash`): as
/// `encode_protein_query_frame_with_seg` (`_timed` with `timed`), the SEG intervals kept in
/// `seg_masks` for the lookup table and the reports, but the letters left as they are (no X, no
/// `aa_seq_nomask`), so the extensions, the Karlin blocks and the composition read the query.
pub fn encode_protein_query_frame_soft_masked(
    seq: &[u8],
    seg: Option<&SegParams>,
    timed: bool,
) -> QueryFrame {
    encode_query_frame::<true>(seq, seg, timed)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1239-1256
// ```c
//     if (filter_options->segOptions) { ... SeqBufferSeg(...) ... }
//     ...
//         BlastSeqLoc *lcase_mask_slp = NULL;
//         if (query_blk->lcase_mask && query_blk->lcase_mask->seqloc_array)
//         {
//             ASSERT(context < query_blk->lcase_mask->total_size);
//             lcase_mask_slp = query_blk->lcase_mask->seqloc_array[context];
//             /* Set location list to NULL to allow safe memory deallocation,
//               ownership transferred to filter_out. */
//             query_blk->lcase_mask->seqloc_array[context] = NULL;
//         }
//         BlastSeqLocAppend(filter_out, lcase_mask_slp);
//         status = BlastSeqLocCombine(filter_out, 0);
// ```
/// The query frame of a BLASTP search with `-lcase_masking`: its SEG intervals (none without
/// `-seg`) and its lower-case runs in `seg_masks`, combined as NCBI combines them. Without `soft`
/// the frame keeps an unmasked copy and its masked letters become X. Cold and out of line (no NCBI
/// counterpart: code placement only), so that the callers keep the default path's code.
#[cold]
#[inline(never)]
pub fn encode_protein_query_frame_lcase(
    seq: &[u8],
    seg: Option<&SegParams>,
    soft: bool,
    timed: bool,
) -> QueryFrame {
    let mut frame = encode_query_frame::<true>(seq, seg, timed);
    // SEG's intervals first, as `BlastSeqLocAppend(filter_out, lcase_mask_slp)`: the sort of
    // `BlastSeqLocCombine` compares only the left ends, and a one-letter run on the first letter of
    // a SEG interval is merged only in this order (the results' masks).
    let mut masks: Vec<(i32, i32)> = frame
        .seg_masks
        .iter()
        .map(|&(start, end)| (start as i32, end as i32 - 1))
        .collect();
    masks.extend(lowercase_runs(seq));
    frame.seg_masks = combine_masks(masks)
        .into_iter()
        .map(|(from, to)| (from as usize, to as usize + 1))
        .collect();
    if !soft && !frame.seg_masks.is_empty() {
        frame.aa_seq_nomask = Some(frame.aa_seq.clone());
        for &(start, end) in &frame.seg_masks {
            // Past the leading sentinel.
            frame.aa_seq[1 + start..1 + end].fill(ncbistdaa::X);
        }
    }
    frame
}

// `MASK_AT_HASH` is a const parameter so that the default (hard) instance and its callers keep the
// code they had before `-soft_masking` (one instance per value).
fn encode_query_frame<const MASK_AT_HASH: bool>(
    seq: &[u8],
    seg: Option<&SegParams>,
    timed: bool,
) -> QueryFrame {
    // Diagnostic (LOSAT_TIMING): the encoding for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_encode = if timed {
        crate::utils::xtiming::clock()
    } else {
        None
    };
    let encoded = encode_protein_sequence(seq);
    let mut aa_seq = encoded.aa_seq.clone();
    let mut aa_seq_nomask = None;
    let mut seg_masks = Vec::new();
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Encode, t_encode);

    if let Some(params) = seg {
        // Diagnostic (LOSAT_TIMING): the query's SEG for the common line (utils/xtiming.rs). No
        // NCBI counterpart.
        let t_seg = if timed {
            crate::utils::xtiming::clock()
        } else {
            None
        };
        let masker = SegMasker::with_params(params);
        let intervals = masker.mask_sequence(&aa_seq[1..aa_seq.len() - 1]);
        if MASK_AT_HASH {
            // `BlastSetUp_MaskQuery` is not called: the intervals only keep words out of the
            // lookup table (`BLAST_ComplementMaskLocations`).
            seg_masks.extend(
                intervals
                    .iter()
                    .map(|interval| (interval.start, interval.end)),
            );
        } else if !intervals.is_empty() {
            aa_seq_nomask = Some(aa_seq.clone());
            for interval in &intervals {
                seg_masks.push((interval.start, interval.end));
                for pos in interval.start..interval.end {
                    aa_seq[pos + 1] = ncbistdaa::X;
                }
            }
        }
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
    // ```c
    //     if (!mask_at_hash) {
    //         BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc,
    //                              program_number);
    //     }
    // ```
    // `-soft_masking true` keeps the SEG intervals of the hard query and its letters.
    #[test]
    fn test_soft_masked_query_keeps_letters_and_seg_intervals() {
        let seq = b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKQQQQQQQQQQQQQQQQQQQQQGWENYCWPANMSWKVLMKYEFSSD";
        let seg = SegParams::default();
        let hard = encode_protein_query_frame_with_seg(seq, Some(&seg));
        let soft = encode_protein_query_frame_soft_masked(seq, Some(&seg), false);
        assert!(!hard.seg_masks.is_empty());
        assert_eq!(soft.seg_masks, hard.seg_masks);
        assert_eq!(soft.aa_seq, encode_protein_sequence(seq).aa_seq);
        assert_eq!(hard.aa_seq_nomask.as_deref(), Some(&soft.aa_seq[..]));
        assert!(soft.aa_seq_nomask.is_none());
        assert_ne!(soft.aa_seq, hard.aa_seq);
        let none = encode_protein_query_frame_soft_masked(seq, None, false);
        assert!(none.seg_masks.is_empty() && none.aa_seq == soft.aa_seq);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1253-1254
    // ```c
    //         BlastSeqLocAppend(filter_out, lcase_mask_slp);
    //         status = BlastSeqLocCombine(filter_out, 0);
    // ```
    // A lower-case run joins the SEG intervals; hard writes X and keeps the letters in the copy,
    // soft keeps the letters.
    #[test]
    fn test_lcase_query_masks_join_seg_intervals() {
        let upper =
            b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKQQQQQQQQQQQQQQQQQQQQQGWENYCWPANMSWKVLMKYEFSSD";
        let mut seq = upper.to_vec();
        seq[2..6].make_ascii_lowercase();
        let seg = SegParams::default();
        let letters = encode_protein_sequence(upper).aa_seq;

        let hard = encode_protein_query_frame_lcase(&seq, Some(&seg), false, false);
        let without = encode_protein_query_frame_with_seg(upper, Some(&seg));
        let mut expected = vec![(2, 6)];
        expected.extend_from_slice(&without.seg_masks);
        assert_eq!(hard.seg_masks, expected);
        assert_eq!(hard.aa_seq_nomask.as_deref(), Some(&letters[..]));
        assert_eq!(hard.aa_seq[3..7], [ncbistdaa::X; 4]);

        let soft = encode_protein_query_frame_lcase(&seq, Some(&seg), true, false);
        assert_eq!(soft.seg_masks, expected);
        assert_eq!(soft.aa_seq, letters);
        assert!(soft.aa_seq_nomask.is_none());

        // No lower case and no SEG: no copy, as without the flag.
        let none = encode_protein_query_frame_lcase(upper, None, false, false);
        assert!(none.seg_masks.is_empty() && none.aa_seq_nomask.is_none());
        assert_eq!(none.aa_seq, letters);
    }

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
