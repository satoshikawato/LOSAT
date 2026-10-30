//! No NCBI counterpart: the preliminary (score-only) greedy gapped alignment of megablast run on
//! a one-byte-per-base copy of the packed ncbi2na subject, with `compressed_subject = false`, so
//! that the uncompressed scans of `BLAST_GreedyAlign` serve the preliminary stage too. A
//! score-only call on the copy runs `x_greedy_fused`, and `x_greedy_row` only where the fused
//! loop declines the call.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4019-4025
//! ```c
//!             status = BLAST_GreedyGappedAlignment(
//!                          query_tmp.sequence, subject->sequence,
//!                          query_tmp.length, subject->length,
//!                          gap_align, score_params,
//!                          init_hsp->offsets.qs_offsets.q_off,
//!                          init_hsp->offsets.qs_offsets.s_off,
//!                          (Boolean) TRUE, FALSE, fence_hit);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2786-2794,2833-2834
//! ```c
//!    if (!compressed_subject) {
//!       s = subject + s_off;
//!       rem = 4; /* Special value to indicate that sequence is already
//!                   uncompressed */
//!    } else {
//!       s = subject + s_off/4;
//!       rem = s_off % 4;
//!    }
//! ...
//!    if (compressed_subject)
//!       rem = 0;
//! ```
//! NCBI runs the preliminary greedy alignment on the compressed subject (`rem` = `s_off % 4`
//! to the right, 0 to the left) and the traceback on the uncompressed one (`rem` = 4). LOSAT
//! gives the preliminary call the uncompressed form too, and has no packed form of the call
//! (it was checked against the packed form when the copy was introduced and deleted at
//! graduation, switch-tidy D2): byte `i` of the copy is `NCBI2NA_UNPACK_BASE(packed[i / 4],
//! 3 - i % 4)`, the base the compressed scan reads. The scan predicate is the same (a query byte
//! `>= 4` never matches, a byte below 4 matches iff it equals the base); the copy holds only
//! 0..3, so it never reports a fence, and the compressed scan has none. To the right the call
//! starts at base `s_off` of the chunk (the chunk starts at a multiple of 4, `next_chunk` in
//! `run.rs`), to the left it reads base `s_off - 1 - i` of the whole chunk slice. The
//! traceback blocks and `reduce_gaps`, which read subject bytes, run only with `do_traceback`,
//! which is false here.

use super::greedy::{greedy_gapped_alignment_internal, GreedyAlignScratch};
use crate::core::blast_encoding::COMPRESSION_RATIO;
use std::sync::OnceLock;

/// The four bases of each packed byte, first base first.
const X_UNPACKED: [[u8; COMPRESSION_RATIO]; 256] = {
    let mut table = [[0u8; COMPRESSION_RATIO]; 256];
    let mut byte = 0usize;
    while byte < 256 {
        let b = byte as u8;
        table[byte] = [(b >> 6) & 3, (b >> 4) & 3, (b >> 2) & 3, b & 3];
        byte += 1;
    }
    table
};

/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_util.h:52-55
/// ```c
/// #define NCBI2NA_MASK 0x03
/// #define NCBI2NA_UNPACK_BASE(x, N) (((x)>>(2*(N))) & NCBI2NA_MASK)
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:365-368
/// ```c
///             while (seq1_index < len1 && seq2_index < len2 &&
///                 seq1[seq1_index] ==
///                 NCBI2NA_UNPACK_BASE(seq2[(seq2_index + rem) / 4],
///                                      3 - (seq2_index + rem) % 4)) {
/// ```
/// Byte `i` of the result is `NCBI2NA_UNPACK_BASE(packed[i / 4], 3 - i % 4)` for the first
/// `len` bases (the remainder count in the last packed byte is not a base and is not copied).
pub(crate) fn x_unpack_ncbi2na(packed: &[u8], len: usize) -> Vec<u8> {
    let bytes = len.div_ceil(COMPRESSION_RATIO);
    assert!(
        packed.len() >= bytes,
        "packed subject shorter than its length"
    );
    let mut unpacked = vec![0u8; bytes * COMPRESSION_RATIO];
    for (bases, &byte) in unpacked
        .chunks_exact_mut(COMPRESSION_RATIO)
        .zip(&packed[..bytes])
    {
        bases.copy_from_slice(&X_UNPACKED[byte as usize]);
    }
    unpacked.truncate(len);
    unpacked
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
/// ```c
///         SBlastSequence compressed_seq =
///             subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
///                                       eNa_strand_plus, eNoSentinels);
///         BlastSeqBlkSetCompressedSequence(subj,
///                                   compressed_seq.data.release());
/// ```
/// The one-byte copy of one subject's packed sequence, made on the first chunk that reaches the
/// gapped stage and shared by every chunk and worker of that subject. Off (`None` from
/// `chunk`) unless this is the preliminary stage of a search that uses the greedy aligner.
/// No NCBI counterpart for the copy itself: NCBI reads the packed bytes.
pub(crate) struct XPrelimSubject<'a> {
    packed_full: &'a [u8],
    len_full: usize,
    enabled: bool,
    unpacked: OnceLock<Vec<u8>>,
}

impl<'a> XPrelimSubject<'a> {
    /// `greedy_prelim`: the preliminary stage of a search that uses the greedy aligner.
    pub(crate) fn new(greedy_prelim: bool, packed_full: &'a [u8], len_full: usize) -> Self {
        Self {
            packed_full,
            len_full,
            enabled: greedy_prelim,
            unpacked: OnceLock::new(),
        }
    }

    /// The copy from base `chunk_offset` (a multiple of 4) to the end of the subject: the bases
    /// of `packed_full[chunk_offset / 4..]`. `None` when the copy is off.
    pub(crate) fn chunk(&self, chunk_offset: usize) -> Option<&[u8]> {
        if !self.enabled {
            return None;
        }
        debug_assert_eq!(chunk_offset % COMPRESSION_RATIO, 0);
        let unpacked = self
            .unpacked
            .get_or_init(|| x_unpack_ncbi2na(self.packed_full, self.len_full));
        unpacked.get(chunk_offset..)
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2762-2794
/// ```c
/// BLAST_GreedyGappedAlignment(const Uint1* query, const Uint1* subject,
///    Int4 query_length, Int4 subject_length, BlastGapAlignStruct* gap_align,
///    const BlastScoringParameters* score_params,
///    Int4 q_off, Int4 s_off, Boolean compressed_subject, Boolean do_traceback,
///    Boolean * fence_hit)
/// ...
///    if (!compressed_subject) {
///       s = subject + s_off;
///       rem = 4; /* Special value to indicate that sequence is already
///                   uncompressed */
/// ```
/// The preliminary greedy alignment (the call at blast_gapalign.c:4019-4025), on the one-byte
/// copy of the chunk (`XPrelimSubject::chunk`) with `compressed_subject = false` (module doc:
/// same result as the packed call). LOSAT has no packed form of this call.
#[allow(clippy::too_many_arguments)]
pub fn x_greedy_gapped_alignment_score_only(
    subject_unpacked: &[u8],
    query: &[u8],
    subject_len: usize,
    q_off: usize,
    s_off: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut GreedyAlignScratch,
) -> Option<(usize, usize, usize, usize, i32, usize, usize)> {
    let core = greedy_gapped_alignment_internal(
        query,
        subject_unpacked,
        subject_len,
        q_off,
        s_off,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop,
        scratch,
        false,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2539-2569
    // (s_BlastGreedyGapAlignStructFill: query_start/stop, subject_start/stop, score)
    let core = core?;
    if core.q_start < 0 || core.s_start < 0 {
        return None;
    }
    Some((
        core.q_start as usize,
        core.q_end as usize,
        core.s_start as usize,
        core.s_end as usize,
        core.score,
        core.q_seed_start as usize,
        core.s_seed_start as usize,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::blast_encoding::{
        encode_subject_ncbi2na_packed, iupacna_to_ncbi4na, resolve_ncbi4na_to_ncbi2na,
    };

    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_util.h:52-55
    /// ```c
    /// #define NCBI2NA_MASK 0x03
    /// #define NCBI2NA_UNPACK_BASE(x, N) (((x)>>(2*(N))) & NCBI2NA_MASK)
    /// ```
    /// The one-byte copy of a packed subject (IUPAC letters resolved as the engine packs them)
    /// holds the resolved base at every position, and `XPrelimSubject::chunk` returns the tail
    /// that starts at a multiple of 4.
    #[test]
    fn x_unpack_ncbi2na_gives_the_resolved_bases() {
        const LETTERS: &[u8] = b"ACGTacgtNRYKMSWBDHVnrykmswbdhv-";
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        for len in [1usize, 3, 4, 5, 63, 64, 65, 1000] {
            let text: Vec<u8> = (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    LETTERS[(state % LETTERS.len() as u64) as usize]
                })
                .collect();
            let packed = encode_subject_ncbi2na_packed(&text);
            let codes = resolve_ncbi4na_to_ncbi2na(
                &text
                    .iter()
                    .map(|&letter| iupacna_to_ncbi4na(letter).unwrap_or(15))
                    .collect::<Vec<u8>>(),
            );
            assert_eq!(x_unpack_ncbi2na(&packed, len), codes, "length {len}");
            let subject = XPrelimSubject::new(true, &packed, len);
            if len >= 4 {
                assert_eq!(subject.chunk(4).expect("enabled"), &codes[4..]);
            }
            assert!(XPrelimSubject::new(false, &packed, len).chunk(0).is_none());
        }
    }
}
