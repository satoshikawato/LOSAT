//! Translation and frame generation for TBLASTX
//!
//! This module handles translation of DNA sequences to amino acid sequences
//! and generation of all six reading frames.
//!
//! NCBI BLAST style: Each translated sequence is wrapped with sentinel bytes
//! at the beginning and end to ensure proper extension termination.
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:438-453

use super::constants::SENTINEL_BYTE;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::aa_char_to_ncbistdaa;
use bio::alphabets::dna;

/// A query frame containing translated sequence and metadata
#[derive(Debug, Clone)]
pub struct QueryFrame {
    /// Frame number: 1..3 for forward, -1..-3 for reverse
    pub frame: i8,
    /// Translated amino acid sequence with NCBI-style sentinel bytes.
    ///
    /// Layout: [SENTINEL_BYTE, aa_0, aa_1, ..., aa_n-1, SENTINEL_BYTE]
    /// - aa_seq[0] = SENTINEL_BYTE (0) - NCBI NULLB
    /// - aa_seq[1..len-1] = actual amino acids (encoded in NCBISTDAA 1-27)
    /// - aa_seq[len-1] = SENTINEL_BYTE (0) - NCBI NULLB
    ///
    /// When accessing amino acid at logical position i, use aa_seq[i + 1].
    /// When reporting coordinates, subtract 1 from the raw position.
    ///
    /// Reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:438-453
    ///   prot_seq[0] = NULLB;  // sentinel at beginning
    ///   prot_seq[index_prot] = NULLB;  // sentinel at end
    pub aa_seq: Vec<u8>,
    /// Unmasked amino acid sequence (NCBI `sequence_nomask` semantics), including sentinels.
    ///
    /// NCBI BLAST preserves a copy of the unmasked query sequence buffer, then overwrites
    /// masked residues in the working sequence with X.
    ///
    /// NCBI reference (verbatim):
    ///   query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
    ///   query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
    pub aa_seq_nomask: Option<Vec<u8>>,
    /// Number of actual amino acids (excluding sentinels)
    pub aa_len: usize,
    /// Original DNA sequence length
    pub orig_len: usize,
    /// SEG-masked AA intervals (logical positions, frame-specific)
    /// Seeds should not be generated from these positions
    pub seg_masks: Vec<(usize, usize)>,
}

/// Convert a codon to amino acid index in NCBISTDAA encoding (0-27)
/// Reference: NCBI BLAST uses NCBISTDAA for all internal sequence representation
fn codon_to_aa_ncbistdaa(codon: &[u8], table: &GeneticCode) -> u8 {
    let aa = table.get(codon);
    aa_char_to_ncbistdaa(aa)
}

/// Translate a DNA sequence to an amino acid sequence with NCBI-style sentinels.
///
/// NCBI BLAST places NULLB sentinel bytes at the beginning and end of each
/// translated sequence to ensure proper extension termination.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:438-453
///   prot_seq[0] = NULLB;  // First character is sentinel
///   prot_seq[index_prot] = NULLB;  // Last character is sentinel
///
/// Returns (aa_seq_with_sentinels, actual_aa_count)
fn translate_sequence(seq: &[u8], table: &GeneticCode) -> (Vec<u8>, usize) {
    let len_aa = seq.len() / 3;
    // +2 for sentinel bytes at beginning and end
    let mut aa_seq = Vec::with_capacity(len_aa + 2);

    // NCBI BLAST style: sentinel at the beginning
    aa_seq.push(SENTINEL_BYTE);

    for chunk in seq.chunks(3) {
        if chunk.len() == 3 {
            aa_seq.push(codon_to_aa_ncbistdaa(chunk, table));
        }
    }

    // NCBI BLAST style: sentinel at the end
    aa_seq.push(SENTINEL_BYTE);

    (aa_seq, len_aa)
}

/// Generate all six reading frames (3 forward + 3 reverse) for a DNA sequence
///
/// # Arguments
/// * `seq` - DNA sequence to translate
/// * `table` - Genetic code table to use for translation
///
/// # Returns
/// Vector of QueryFrame structures, one for each valid reading frame.
/// Each frame's aa_seq contains sentinel bytes at positions 0 and len-1.
/// Amino acids are encoded in NCBISTDAA (0-27).
pub fn generate_frames(seq: &[u8], table: &GeneticCode) -> Vec<QueryFrame> {
    let mut frames = Vec::with_capacity(6);
    let rev_seq = dna::revcomp(seq);
    let seq_len = seq.len();

    // Forward (Frames 1,2,3)
    for i in 0..3 {
        if i + 3 <= seq_len {
            let (aa_seq, aa_len) = translate_sequence(&seq[i..], table);
            frames.push(QueryFrame {
                frame: (i as i8) + 1,
                aa_seq,
                aa_seq_nomask: None,
                aa_len,
                orig_len: seq_len,
                seg_masks: Vec::new(),
            });
        }
    }
    // Reverse (Frames -1,-2,-3)
    for i in 0..3 {
        if i + 3 <= rev_seq.len() {
            let (aa_seq, aa_len) = translate_sequence(&rev_seq[i..], table);
            frames.push(QueryFrame {
                frame: -((i as i8) + 1),
                aa_seq,
                aa_seq_nomask: None,
                aa_len,
                orig_len: seq_len,
                seg_masks: Vec::new(),
            });
        }
    }
    frames
}

/// The bases a window of `x_generate_frames_in_windows` holds: one window of a subject's letters
/// and two letters of the window before it (3 x 2^16 bases, 192 KiB).
const X_FRAME_WINDOW: usize = 3 << 16;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1296-1308
// ```c
//    for (context = 0; context < num_frames; ++context) {
//       int frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//       retval->translations[context] = (Uint1*) malloc(...);
//       BLAST_GetTranslation(subject_blk->sequence_start, nucl_seq_rev,
//          subject_blk->length, frame, retval->translations[context], gen_code_string);
//    }
// ```
// EXPERIMENT (LOSAT_X_TBXFRAMES): `generate_frames` of a subject whose letters are read in
// windows, without the subject's letters or their reverse complement in memory at once: the
// six frames are filled in one pass over the letters, each codon start q (0 <= q <= len - 3)
// giving the residue of forward frame q % 3 at index q / 3 (the codon letters q..q+3) and of
// reverse frame r % 3 at index r / 3 with r = len - 3 - q (the complements of letters q+2, q+1,
// q, in that order: `dna::revcomp` is `rev().map(complement)`). Every residue of the six frames
// of `generate_frames` is the translation of the same three letters by the same
// `codon_to_aa_ncbistdaa`, so the frames are the same: the frames that exist (`i + 3 <= len`),
// in the same order, with the same lengths, sentinels and residues.
/// The frames `generate_frames` makes of the `len` letters that `fill` writes: `fill(from, out)`
/// writes the letters `[from, from + out.len())`; it is called for consecutive windows that cover
/// `[0, len)` once, in order.
pub(crate) fn x_generate_frames_in_windows(
    len: usize,
    fill: &mut dyn FnMut(usize, &mut [u8]),
    table: &GeneticCode,
) -> Vec<QueryFrame> {
    x_frames_in_windows(len, X_FRAME_WINDOW, fill, table)
}

fn x_frames_in_windows(
    len: usize,
    window: usize,
    fill: &mut dyn FnMut(usize, &mut [u8]),
    table: &GeneticCode,
) -> Vec<QueryFrame> {
    let frame_of = |frame: i8, i: usize| {
        let aa_len = (len - i) / 3;
        QueryFrame {
            frame,
            aa_seq: vec![SENTINEL_BYTE; aa_len + 2],
            aa_seq_nomask: None,
            aa_len,
            orig_len: len,
            seg_masks: Vec::new(),
        }
    };
    let mut forward: Vec<QueryFrame> = (0..3)
        .filter(|i| i + 3 <= len)
        .map(|i| frame_of(i as i8 + 1, i))
        .collect();
    let mut reverse: Vec<QueryFrame> = (0..3)
        .filter(|i| i + 3 <= len)
        .map(|i| frame_of(-(i as i8 + 1), i))
        .collect();
    // `letters` holds the letters [start, end): the codon starts not yet translated (at most two
    // letters of the last window) and the window just filled.
    let mut letters: Vec<u8> = Vec::with_capacity(window.min(len) + 2);
    let mut start = 0usize;
    let mut end = 0usize;
    while end < len {
        let to = (end + window).min(len);
        let kept = letters.len();
        letters.resize(kept + (to - end), 0);
        fill(end, &mut letters[kept..]);
        end = to;
        // The codon starts whose three letters are in [start, end).
        let codons_end = end.saturating_sub(2).max(start);
        for q in start..codons_end {
            let at = q - start;
            let (b0, b1, b2) = (letters[at], letters[at + 1], letters[at + 2]);
            forward[q % 3].aa_seq[q / 3 + 1] = codon_to_aa_ncbistdaa(&[b0, b1, b2], table);
            let r = len - 3 - q;
            reverse[r % 3].aa_seq[r / 3 + 1] = codon_to_aa_ncbistdaa(
                &[
                    dna::complement(b2),
                    dna::complement(b1),
                    dna::complement(b0),
                ],
                table,
            );
        }
        letters.drain(..codons_end - start);
        start = codons_end;
    }
    forward.append(&mut reverse);
    forward
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `x_frames_in_windows` against `generate_frames`: every field of every frame.
    fn assert_same_frames(seq: &[u8], window: usize, code: &GeneticCode) {
        let expected = generate_frames(seq, code);
        let mut fills: Vec<(usize, usize)> = Vec::new();
        let got = x_frames_in_windows(
            seq.len(),
            window,
            &mut |from, out| {
                fills.push((from, out.len()));
                out.copy_from_slice(&seq[from..from + out.len()]);
            },
            code,
        );
        assert_eq!(
            got.len(),
            expected.len(),
            "len {} window {window}",
            seq.len()
        );
        for (g, e) in got.iter().zip(&expected) {
            assert_eq!(g.frame, e.frame);
            assert_eq!(
                g.aa_seq,
                e.aa_seq,
                "frame {} len {} window {window}",
                e.frame,
                seq.len()
            );
            assert!(g.aa_seq_nomask.is_none() && e.aa_seq_nomask.is_none());
            assert_eq!(g.aa_len, e.aa_len);
            assert_eq!(g.orig_len, e.orig_len);
            assert!(g.seg_masks.is_empty() && e.seg_masks.is_empty());
        }
        // The windows cover the letters once, in order.
        let mut next = 0;
        for (from, n) in fills {
            assert_eq!(from, next);
            assert!(n > 0);
            next += n;
        }
        assert_eq!(next, seq.len());
    }

    #[test]
    fn frames_in_windows_are_the_frames_of_the_letters() {
        // Letters of both cases, every IUPAC code, U, a gap and bytes that are not letters.
        let alphabet = b"ACGTACGTACGTacgtNnRYKMSWBDHVrykmswbdhvUu-*X.";
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut letter = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            alphabet[(state % alphabet.len() as u64) as usize]
        };
        let codes = [1u8, 2, 4, 11];
        for len in (0..40).chain([95, 96, 97, 98, 99, 100, 101, 1000, 3001]) {
            let seq: Vec<u8> = (0..len).map(|_| letter()).collect();
            for window in [1, 2, 3, 4, 5, 7, 32, 33, 34, 100, X_FRAME_WINDOW] {
                for &id in &codes {
                    assert_same_frames(&seq, window, &GeneticCode::from_id(id));
                }
            }
        }
        let long: Vec<u8> = (0..2 * X_FRAME_WINDOW + 5).map(|_| letter()).collect();
        assert_same_frames(&long, X_FRAME_WINDOW, &GeneticCode::from_id(1));
    }

    #[test]
    fn test_generate_frames() {
        let code = GeneticCode::from_id(1);
        let seq = b"ATGCGATCGATCGATCG";
        let frames = generate_frames(seq, &code);

        // Should generate up to 6 frames
        assert!(frames.len() <= 6);

        // Check that frames have valid frame numbers and sentinel structure
        for frame in &frames {
            assert!(frame.frame >= -3 && frame.frame <= 3 && frame.frame != 0);
            assert_eq!(frame.orig_len, seq.len());
            // Check sentinels
            assert_eq!(frame.aa_seq[0], SENTINEL_BYTE);
            assert_eq!(frame.aa_seq[frame.aa_seq.len() - 1], SENTINEL_BYTE);
            // Check aa_len matches actual amino acids
            assert_eq!(frame.aa_len, frame.aa_seq.len() - 2);
        }
    }

    #[test]
    fn test_sentinel_positions() {
        let code = GeneticCode::from_id(1);
        let seq = b"ATGATGATGATG"; // 4 codons = 4 amino acids
        let frames = generate_frames(seq, &code);

        let frame = &frames[0]; // Frame 1
                                // aa_seq should be: [SENTINEL(0), aa0, aa1, aa2, aa3, SENTINEL(0)]
        assert_eq!(frame.aa_seq.len(), 6); // 4 AAs + 2 sentinels
        assert_eq!(frame.aa_len, 4);
        assert_eq!(frame.aa_seq[0], SENTINEL_BYTE);
        assert_eq!(frame.aa_seq[5], SENTINEL_BYTE);
        // Middle positions should be valid amino acids (not sentinels)
        // Note: NCBISTDAA encodes actual AAs as 1-27, so they are never 0
        for i in 1..5 {
            assert_ne!(frame.aa_seq[i], SENTINEL_BYTE);
        }
    }

    #[test]
    fn test_ncbistdaa_encoding() {
        use crate::utils::matrix::ncbistdaa;

        let code = GeneticCode::from_id(1);
        // ATG = Met (M), NCBISTDAA M = 12
        let seq = b"ATG";
        let frames = generate_frames(seq, &code);
        let frame = &frames[0];
        assert_eq!(frame.aa_seq[1], ncbistdaa::M);

        // TAA = Stop (*), NCBISTDAA * = 25
        let seq_stop = b"TAA";
        let frames_stop = generate_frames(seq_stop, &code);
        let frame_stop = &frames_stop[0];
        assert_eq!(frame_stop.aa_seq[1], ncbistdaa::STOP);
    }
}
