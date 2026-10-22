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

/// The residue of codon `offset` of `frame` (+1..+3, -1..-3) of `sequence` as the search
/// translates it: the residue `generate_frames` puts at index `offset + 1` of that frame (the
/// codon from `3 * offset + |frame| - 1` of the sequence, or of its reverse complement
/// `dna::revcomp`, translated by `codon_to_aa_ncbistdaa`). The report reads single codons
/// with it (the engine counts of the Seq-aligns, `report::engine_counts`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1299-1308
/// ```c
///                /* First produce the reverse strand of the nucleotide sequence */
///               GetReverseNuclSequence(subject_blk->sequence_start, subject_blk->length,
///                     &nucl_seq_rev);
///
///               for (context = 0; context < num_frames; ++context) {
///                  int frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
///                  retval->translations[context] = (Uint1*) malloc((2+subject_blk->length/3)*sizeof(Uint1));
///                  BLAST_GetTranslation(subject_blk->sequence_start, nucl_seq_rev,
///                     subject_blk->length, frame, retval->translations[context], gen_code_string);
///               }
/// ```
pub(crate) fn frame_residue(sequence: &[u8], frame: i8, offset: usize, table: &GeneticCode) -> u8 {
    let first = 3 * offset + frame.unsigned_abs() as usize - 1;
    if frame > 0 {
        codon_to_aa_ncbistdaa(&sequence[first..first + 3], table)
    } else {
        // `dna::revcomp(sequence)[j]` is `dna::complement(sequence[len - 1 - j])`.
        let last = sequence.len() - 1 - first;
        codon_to_aa_ncbistdaa(
            &[
                dna::complement(sequence[last]),
                dna::complement(sequence[last - 1]),
                dna::complement(sequence[last - 2]),
            ],
            table,
        )
    }
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1216-1218
// ```c
//                length = BLAST_GetTranslation(nucl_seq, nucl_seq_rev,
//                        nucl_length, hsp->subject.frame, target_t->translations[context],
//                        target_t->gen_code_string);
// ```
// NCBI translates the one frame of the HSP. The frame here is the one `generate_frames` makes:
// the same `translate_sequence` of the same letters (`seq[i..]`, or `revcomp(seq)[i..]`), present
// under the same condition (`i + 3 <= len`).
/// The frame `frame` (+1..+3, -1..-3) of `generate_frames(seq, table)`, as (residues with the
/// sentinels, length), or `None` where `generate_frames` makes no such frame. Only that frame is
/// translated (and for a reverse frame the reverse complement made).
pub(crate) fn translate_frame(
    seq: &[u8],
    frame: i8,
    table: &GeneticCode,
) -> Option<(Vec<u8>, usize)> {
    let i = usize::from(frame.unsigned_abs()).checked_sub(1)?;
    if i >= 3 || i + 3 > seq.len() {
        return None;
    }
    Some(if frame > 0 {
        translate_sequence(&seq[i..], table)
    } else {
        translate_sequence(&dna::revcomp(seq)[i..], table)
    })
}

/// The bases a window of `x_generate_frames_in_windows` holds: one window of a subject's letters
/// and two letters of the window before it (3 x 2^16 bases, 192 KiB).
const X_FRAME_WINDOW: usize = 3 << 16;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1299-1308
// ```c
//                /* First produce the reverse strand of the nucleotide sequence */
//               GetReverseNuclSequence(subject_blk->sequence_start, subject_blk->length,
//                     &nucl_seq_rev);
//
//               for (context = 0; context < num_frames; ++context) {
//                  int frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//                  retval->translations[context] = (Uint1*) malloc((2+subject_blk->length/3)*sizeof(Uint1));
//                  BLAST_GetTranslation(subject_blk->sequence_start, nucl_seq_rev,
//                     subject_blk->length, frame, retval->translations[context], gen_code_string);
//               }
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
    let residues = x_codon_residues(table);
    let masks = x_letter_masks();
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
        if forward.len() == 3 {
            // Codon starts q = q0, q0 + 3, ... of forward frame f (q % 3 = f, index q / 3, ascending)
            // are those of one reverse frame ((len - 3 - q) % 3, index (len - 3 - q) / 3, descending).
            for f in 0..3 {
                let q0 = start + (f + 3 - start % 3) % 3;
                if q0 >= codons_end {
                    continue;
                }
                let count = (codons_end - q0).div_ceil(3);
                let k0 = q0 / 3;
                let r0 = len - 3 - q0;
                let (ir, kr0) = (r0 % 3, r0 / 3);
                let fw = &mut forward[f].aa_seq[1 + k0..1 + k0 + count];
                let rv = &mut reverse[ir].aa_seq[2 + kr0 - count..2 + kr0];
                let codons = letters[q0 - start..end - start].chunks_exact(3);
                for ((codon, out_f), out_r) in codons.zip(fw.iter_mut()).zip(rv.iter_mut().rev()) {
                    let (b0, b1, b2) = (codon[0] as usize, codon[1] as usize, codon[2] as usize);
                    *out_f = residues[usize::from(masks.forward[b0]) << 8
                        | usize::from(masks.forward[b1]) << 4
                        | usize::from(masks.forward[b2])];
                    *out_r = residues[usize::from(masks.reverse[b2]) << 8
                        | usize::from(masks.reverse[b1]) << 4
                        | usize::from(masks.reverse[b0])];
                }
            }
        } else {
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
        }
        letters.drain(..codons_end - start);
        start = codons_end;
    }
    forward.append(&mut reverse);
    forward
}

/// The ncbi4na mask (`base_mask`) of every byte as a codon letter, and of its complement
/// (`dna::complement`, the map of `dna::revcomp`).
struct XLetterMasks {
    forward: [u8; 256],
    reverse: [u8; 256],
}

fn x_letter_masks() -> &'static XLetterMasks {
    static MASKS: std::sync::OnceLock<XLetterMasks> = std::sync::OnceLock::new();
    MASKS.get_or_init(|| XLetterMasks {
        forward: std::array::from_fn(|byte| crate::utils::genetic_code::base_mask(byte as u8)),
        reverse: std::array::from_fn(|byte| {
            crate::utils::genetic_code::base_mask(dna::complement(byte as u8))
        }),
    })
}

/// A letter of each ncbi4na mask: `base_mask(X_MASK_LETTERS[m]) == m`.
const X_MASK_LETTERS: [u8; 16] = *b"-ACMGRSVTWYHKDBN";

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:369-424 (`s_CodonToAA`; the
// quote is at `GeneticCode::get`)
// ```c
// static Uint1 mapping[4] = { 8, 2, 1, 4 };
// for (i = 0; i < 4; i++) if (codon[0] & mapping[i])
// ```
// `GeneticCode::get` reads a codon only through the masks of its three letters (`base_mask`;
// `LOSAT_X_CODONFAST` gives the loop's residue for single-bit masks), so the residue of a codon
// is the residue of any codon with the same masks: here, of the letters `X_MASK_LETTERS`.
/// The NCBIstdaa residue `codon_to_aa_ncbistdaa` gives every codon with masks (m0, m1, m2), at
/// index m0 << 8 | m1 << 4 | m2; one table per thread for the last genetic code asked for.
fn x_codon_residues(table: &GeneticCode) -> std::rc::Rc<[u8; 4096]> {
    thread_local! {
        static LAST: std::cell::RefCell<Option<([u8; 64], std::rc::Rc<[u8; 4096]>)>> =
            const { std::cell::RefCell::new(None) };
    }
    LAST.with(|last| {
        let mut last = last.borrow_mut();
        if let Some((code, residues)) = last.as_ref() {
            if *code == table.table {
                return std::rc::Rc::clone(residues);
            }
        }
        let residues = std::rc::Rc::new(std::array::from_fn(|index| {
            let codon = [
                X_MASK_LETTERS[index >> 8],
                X_MASK_LETTERS[(index >> 4) & 15],
                X_MASK_LETTERS[index & 15],
            ];
            codon_to_aa_ncbistdaa(&codon, table)
        }));
        *last = Some((table.table, std::rc::Rc::clone(&residues)));
        residues
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `translate_frame` against `generate_frames`: each frame's residues and length, and no
    /// frame where `generate_frames` makes none (also frames 0 and 4), for every length up to 40 and a long sequence,
    /// with every byte the readers give (IUPAC letters, `-`, both cases), under two codes.
    #[test]
    fn translate_frame_is_the_frame_of_generate_frames() {
        let letters = b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn";
        let mut state = 7u64;
        let mut sequence = |length: usize| -> Vec<u8> {
            (0..length)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    let draw = (state >> 33) as usize;
                    // mostly unambiguous bases, so stop codons and ordinary residues occur
                    if draw % 4 == 0 {
                        letters[draw / 4 % letters.len()]
                    } else {
                        b"ACGTacgt"[draw / 4 % 8]
                    }
                })
                .collect()
        };
        let mut compared = 0;
        for id in [1, 4] {
            let code = GeneticCode::try_from_id(id).unwrap();
            for length in (0..=40).chain([3_001]) {
                let seq = sequence(length);
                let frames = generate_frames(&seq, &code);
                for frame in [1, 2, 3, -1, -2, -3, 0, 4, -4] {
                    let expected = frames
                        .iter()
                        .find(|candidate| candidate.frame == frame)
                        .map(|candidate| (candidate.aa_seq.clone(), candidate.aa_len));
                    assert_eq!(
                        translate_frame(&seq, frame, &code),
                        expected,
                        "{length} {frame}"
                    );
                    compared += usize::from(expected.is_some());
                }
            }
        }
        // lengths 3 and 4 have two and four frames, 5 to 40 and 3,001 all six
        assert_eq!(compared, 2 * (2 + 4 + 6 * 36 + 6));
    }

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
    fn the_mask_letters_have_their_masks() {
        for (mask, &letter) in X_MASK_LETTERS.iter().enumerate() {
            assert_eq!(crate::utils::genetic_code::base_mask(letter) as usize, mask);
        }
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
