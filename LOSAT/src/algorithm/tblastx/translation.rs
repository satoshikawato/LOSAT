//! Translation and frame generation for TBLASTX
//!
//! This module handles translation of DNA sequences to amino acid sequences
//! and generation of all six reading frames.
//!
//! NCBI BLAST style: Each translated sequence is wrapped with sentinel bytes
//! at the beginning and end to ensure proper extension termination.
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:438-453

use super::constants::SENTINEL_BYTE;
use crate::utils::genetic_code::{complement, revcomp, GeneticCode};
use crate::utils::matrix::aa_char_to_ncbistdaa;

/// A query frame containing translated sequence and metadata
#[derive(Debug, Clone)]
pub struct QueryFrame {
    /// Frame number: 1..3 for forward, -1..-3 for reverse
    pub frame: i8,
    /// Translated amino acid sequence with NCBI-style sentinel bytes.
    ///
    /// Layout: [SENTINEL_BYTE, aa_0, aa_1, ..., aa_n-1, SENTINEL_BYTE]
    /// - `aa_seq[0]` = SENTINEL_BYTE (0) - NCBI NULLB
    /// - `aa_seq[1..len-1]` = actual amino acids (encoded in NCBISTDAA 1-27)
    /// - `aa_seq[len-1]` = SENTINEL_BYTE (0) - NCBI NULLB
    ///
    /// When accessing amino acid at logical position i, use `aa_seq[i + 1]`.
    /// When reporting coordinates, subtract 1 from the raw position.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:438-453
    /// The first and the last letter of each translated frame are the `NULLB` sentinel.
    pub aa_seq: Vec<u8>,
    /// Unmasked amino acid sequence (NCBI `sequence_nomask` semantics), including sentinels.
    ///
    /// NCBI BLAST preserves a copy of the unmasked query sequence buffer, then overwrites
    /// masked residues in the working sequence with X.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1381-1382
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:438-453
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
/// `revcomp`, translated by `codon_to_aa_ncbistdaa`). The report reads single codons
/// with it (the engine counts of the Seq-aligns, `report::engine_counts`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1299-1308
pub(crate) fn frame_residue(sequence: &[u8], frame: i8, offset: usize, table: &GeneticCode) -> u8 {
    let first = 3 * offset + frame.unsigned_abs() as usize - 1;
    if frame > 0 {
        codon_to_aa_ncbistdaa(&sequence[first..first + 3], table)
    } else {
        // `revcomp(sequence)[j]` is `complement(sequence[len - 1 - j])`.
        let last = sequence.len() - 1 - first;
        codon_to_aa_ncbistdaa(
            &[
                complement(sequence[last]),
                complement(sequence[last - 1]),
                complement(sequence[last - 2]),
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
    let rev_seq = revcomp(seq);
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
// NCBI translates the one frame of the HSP. The frame here is the one `generate_frames` makes,
// present under the same condition (`i + 3 <= len`): `translate_sequence` of `seq[i..]` or of
// `revcomp(seq)[i..]`, whose codon k is the letters `seq[i + 3k..i + 3k + 3]`, or the complements
// of `seq[len - i - 3k - 1]`, `seq[len - i - 3k - 2]` and `seq[len - i - 3k - 3]` (`revcomp` is
// `rev().map(complement)`). `x_codons_into` reads those letters in place and gives each codon the
// residue `codon_to_aa_ncbistdaa` gives it.
/// The number of residues of frame `frame` (+1..+3, -1..-3) of `generate_frames` of `len` letters,
/// or `None` where `generate_frames` makes no such frame.
pub(crate) fn frame_length(len: usize, frame: i8) -> Option<usize> {
    let i = usize::from(frame.unsigned_abs()).checked_sub(1)?;
    (i < 3 && i + 3 <= len).then(|| (len - i) / 3)
}

/// The residues of frame `frame` (+1..+3, -1..-3) of `generate_frames(seq, table)` without its
/// sentinels, appended to `out`, and their number; `None` (and `out` unchanged) where
/// `generate_frames` makes no such frame. `residues` is the table of the code
/// (`x_codon_residues(table)`). Only that frame is translated, and no reverse complement is made.
pub(crate) fn translate_frame_into(
    seq: &[u8],
    frame: i8,
    residues: &[u8; 4096],
    out: &mut Vec<u8>,
) -> Option<usize> {
    let count = frame_length(seq.len(), frame)?;
    let i = usize::from(frame.unsigned_abs()) - 1;
    if frame > 0 {
        x_codons_into(&seq[i..i + 3 * count], false, residues, out);
    } else {
        x_codons_into(
            &seq[seq.len() - i - 3 * count..seq.len() - i],
            true,
            residues,
            out,
        );
    }
    Some(count)
}

/// Appends to `out` the residue (`codon_to_aa_ncbistdaa`, read from the code's table `residues`,
/// `x_codon_residues`) of every codon of `letters` (a multiple of three letters) read forward, or
/// with `reverse` of every codon of their reverse complement (`revcomp(letters)`) read forward.
pub(crate) fn x_codons_into(
    letters: &[u8],
    reverse: bool,
    residues: &[u8; 4096],
    out: &mut Vec<u8>,
) {
    debug_assert_eq!(letters.len() % 3, 0, "whole codons");
    let masks = x_letter_masks();
    if reverse {
        // codon k of `revcomp(letters)` is the complements of the letters of the k-th triple
        // from the end, last letter first
        let mask = |letter: u8| usize::from(masks.reverse[usize::from(letter)]);
        out.extend(
            letters
                .rchunks_exact(3)
                .map(|codon| residues[mask(codon[2]) << 8 | mask(codon[1]) << 4 | mask(codon[0])]),
        );
    } else {
        let mask = |letter: u8| usize::from(masks.forward[usize::from(letter)]);
        out.extend(
            letters
                .chunks_exact(3)
                .map(|codon| residues[mask(codon[0]) << 8 | mask(codon[1]) << 4 | mask(codon[2])]),
        );
    }
}

/// The bases a window of `x_generate_frames_in_windows` holds: one window of a subject's letters
/// and two letters of the window before it (3 x 2^16 bases, 192 KiB).
const X_FRAME_WINDOW: usize = 3 << 16;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1299-1308
// `x_generate_frames_in_windows` is `generate_frames` of a subject whose letters are read in
// windows, without the subject's letters or their reverse complement in memory at once: the
// six frames are filled in one pass over the letters, each codon start q (0 <= q <= len - 3)
// giving the residue of forward frame q % 3 at index q / 3 (the codon letters q..q+3) and of
// reverse frame r % 3 at index r / 3 with r = len - 3 - q (the complements of letters q+2, q+1,
// q, in that order: `revcomp` is `rev().map(complement)`). Every residue of the six frames
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
                reverse[r % 3].aa_seq[r / 3 + 1] =
                    codon_to_aa_ncbistdaa(&[complement(b2), complement(b1), complement(b0)], table);
            }
        }
        letters.drain(..codons_end - start);
        start = codons_end;
    }
    forward.append(&mut reverse);
    forward
}

/// The ncbi4na mask (`base_mask`) of every byte as a codon letter, and of its complement
/// (`complement`, the map of `revcomp`).
struct XLetterMasks {
    forward: [u8; 256],
    reverse: [u8; 256],
}

fn x_letter_masks() -> &'static XLetterMasks {
    static MASKS: std::sync::OnceLock<XLetterMasks> = std::sync::OnceLock::new();
    MASKS.get_or_init(|| XLetterMasks {
        forward: std::array::from_fn(|byte| crate::utils::genetic_code::base_mask(byte as u8)),
        reverse: std::array::from_fn(|byte| {
            crate::utils::genetic_code::base_mask(complement(byte as u8))
        }),
    })
}

/// A letter of each ncbi4na mask: `base_mask(X_MASK_LETTERS[m]) == m`.
const X_MASK_LETTERS: [u8; 16] = *b"-ACMGRSVTWYHKDBN";

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:369-425 (`s_CodonToAA`, which
// `GeneticCode::get` ports)
// `GeneticCode::get` reads a codon only through the masks of its three letters (`base_mask`;
// its table read gives the loop's residue for single-bit masks), so the residue of a codon
// is the residue of any codon with the same masks: here, of the letters `X_MASK_LETTERS`.
/// The NCBIstdaa residue `codon_to_aa_ncbistdaa` gives every codon with masks (m0, m1, m2), at
/// index m0 << 8 | m1 << 4 | m2; one table per thread for the last genetic code asked for.
pub(crate) fn x_codon_residues(table: &GeneticCode) -> std::rc::Rc<[u8; 4096]> {
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
#[path = "translation_inline_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "translation_tests.rs"]
mod translation_tests;
