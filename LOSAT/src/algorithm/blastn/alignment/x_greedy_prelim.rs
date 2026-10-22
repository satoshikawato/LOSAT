//! No NCBI counterpart: the preliminary (score-only) greedy gapped alignment of megablast run on
//! a one-byte-per-base copy of the packed ncbi2na subject, with `compressed_subject = false`, so
//! that the uncompressed scans of `BLAST_GreedyAlign` serve the preliminary stage too. Those are
//! the eight-bases-per-step rows (`x_greedy_row`).
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
//! to the right, 0 to the left) and the traceback on the uncompressed one (`rem` = 4). This
//! module gives the preliminary call the same values in the uncompressed form: byte `i` of the
//! copy is `NCBI2NA_UNPACK_BASE(packed[i / 4], 3 - i % 4)`, the base the compressed scan reads.
//!
//! Why the result is the same (each step checked against the code of `greedy.rs`; `G:` lines
//! are lines of `greedy.rs` at origin/main d6105920):
//!
//! 1. Scan. `find_first_mismatch_greedy` (G:2250-2328, greedy_align.c:313-375) on the compressed
//!    subject continues while `seq1[i] == NCBI2NA_UNPACK_BASE(...)`, a value in 0..3, so a query
//!    byte `seq1 >= 4` (blastna ambiguity or gap) never matches and a byte below 4 matches iff it
//!    equals the base. The uncompressed scan (`rem == 4`) continues while
//!    `seq1[i] < 4 && seq1[i] == seq2[i]` with `seq2[i]` the same base: the same predicate, so
//!    the same run length. The fast rows test `(w1 ^ w2) | (w1 & 0xFC..FC)` per byte (G:2588-2644
//!    and pass 2 at G:2884-2936), the same predicate again.
//! 2. Fence. The uncompressed scan reports a fence when the stop byte of `seq2` is
//!    `FENCE_SENTRY` (201) (G:2275, G:2309; fast rows G:2728 and pass 2). The copy holds only
//!    0..3, so it never does; the compressed scan has no fence test. `fence_hit` stays false
//!    in both.
//! 3. Positions and bounds. To the right the compressed scan reads byte `(i + rem) / 4` of the
//!    slice that starts at byte `s_off / 4`, slot `(i + rem) % 4`, i.e. base `s_off + i` of the
//!    chunk (the chunk starts at a multiple of 4, `next_chunk` in `run.rs`); the uncompressed
//!    scan reads byte `i` of the slice that starts at `s_off`, the same base. To the left both
//!    read base `s_off - 1 - i` of the whole chunk slice (`rem` = 0, `len2 = s_off`). `len1` and
//!    `len2` (`q_avail`, `s_avail`, `q_off`, `s_off`) do not depend on the form. The early return
//!    `subject_offset >= subject.len()` (G:4243) is false in both forms when
//!    `s_off < subject_len`: the packed slice has `len_full / 4 + 1 - offset / 4` bytes and
//!    `s_off / 4 <= (len_full - 1) / 4 - offset / 4`; the copy has `len_full - offset >=
//!    subject_len` bytes. Otherwise `s_avail <= 0` and both forms return `None` (G:4250) before
//!    they touch the scratch. The fast-row dispatch (G:3397-3403) needs `len2 <= seq2.len()`:
//!    to the right the slice has `len_full - offset - s_off >= s_avail` bytes, to the left it is
//!    the whole chunk slice and `len2 = s_off`.
//! 4. Every other use of `rem` / `compressed_subject`: the `rem` and slice choice (G:4229-4249,
//!    item 3), the parameters and the pass-through to `blast_greedy_align` (G:3197, G:3735,
//!    G:3769), the four scans (non-affine d = 0 G:3275 and per diagonal G:3492, affine d = 0
//!    G:3835 and per diagonal G:3969, item 1) each followed by its fence test (item 2), and the
//!    fast-row dispatch `rem == 4` (G:3398); the fast rows give the rows of the ported loop on
//!    `rem == 4` (the test of this module compares the two forms). The traceback blocks and
//!    `reduce_gaps`, which read subject bytes, run only with `do_traceback`, which is false here.
//! 5. Scratch. `max_dist` starts at `max(scratch.max_dist, initial)` with `initial` from
//!    `subject_len` (the same in both forms); a larger `max_dist` only saves a retry.

use super::greedy::{
    greedy_gapped_alignment_internal, greedy_gapped_alignment_score_only, GreedyAlignScratch,
};
use crate::core::blast_encoding::COMPRESSION_RATIO;
use std::sync::OnceLock;

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
    let mut unpacked = Vec::with_capacity(bytes * COMPRESSION_RATIO);
    for &byte in &packed[..bytes] {
        unpacked.extend_from_slice(&[(byte >> 6) & 3, (byte >> 4) & 3, (byte >> 2) & 3, byte & 3]);
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
/// The preliminary greedy alignment (the call at blast_gapalign.c:4019-4025). With the one-byte
/// copy of the chunk (`XPrelimSubject::chunk`) the alignment runs with `compressed_subject =
/// false` on the copy (module doc: same result). With `x_subject = None` (the search has no copy)
/// this is `greedy_gapped_alignment_score_only` on the packed subject, as ported.
#[allow(clippy::too_many_arguments)]
pub fn x_greedy_gapped_alignment_score_only(
    x_subject: Option<&[u8]>,
    query: &[u8],
    subject_packed: &[u8],
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
    let Some(subject_unpacked) = x_subject else {
        return greedy_gapped_alignment_score_only(
            query,
            subject_packed,
            subject_len,
            q_off,
            s_off,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
        );
    };
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
        false,
        scratch,
        false,
    );
    // The tail of `greedy_gapped_alignment_score_only`, unchanged.
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

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:328-370
    /// ```c
    ///         if (rem == 4) {
    ///             while (seq1_index < len1 && seq2_index < len2 &&
    ///                    seq1[len1-1 - seq1_index] < 4 &&
    ///                    seq1[len1-1 - seq1_index] == seq2[len2-1 - seq2_index]) {
    /// ...
    ///             while (seq1_index < len1 && seq2_index < len2 &&
    ///                 seq1[seq1_index] ==
    ///                 NCBI2NA_UNPACK_BASE(seq2[(seq2_index + rem) / 4],
    ///                                      3 - (seq2_index + rem) % 4)) {
    /// ```
    /// On random pairs of related sequences (subjects in IUPAC text with ambiguous letters,
    /// packed as the engine packs them; queries with blastna ambiguity codes in some cases),
    /// seeds at every phase `s_off % 4` and chunks that start at a multiple of 4 after the
    /// subject start and whose slices run past the chunk length: the score-only result on the
    /// one-byte copy (`compressed_subject = false`) equals the result of the packed call as
    /// ported, on new and on used scratch.
    #[test]
    fn x_greedy_prelim_unpacked_matches_packed_on_random_alignments() {
        struct Rng(u64);
        impl Rng {
            fn next(&mut self) -> u64 {
                self.0 ^= self.0 << 13;
                self.0 ^= self.0 >> 7;
                self.0 ^= self.0 << 17;
                self.0
            }
            fn below(&mut self, n: u64) -> u64 {
                self.next() % n
            }
        }
        const BASES: &[u8] = b"ACGTacgt";
        const AMBIGUOUS: &[u8] = b"NRYKMSWBDHVnrykmswbdhv-";

        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(600);
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let mut used = GreedyAlignScratch::new();
        let mut calls = 0usize;
        let mut aligned = 0usize;
        let mut phases = [0usize; 4];
        let mut long_runs = 0usize;
        for case in 0..cases {
            let core_len = 20 + rng.below(if case % 10 == 0 { 6000 } else { 600 }) as usize;
            let substitution = [0u64, 0, 2, 5, 20, 60, 150][rng.below(7) as usize]; // per mille
            let indel = [0u64, 0, 2, 10, 30][rng.below(5) as usize]; // per mille
            let subject_ambiguity = [0u64, 5, 30][rng.below(3) as usize]; // per mille
            let query_ambiguity = if case % 7 == 0 { 8 } else { 0 }; // per mille

            // subject text: flank, core, flank; ambiguous letters resolved by the packing
            let prefix_len = 4 + rng.below(80) as usize;
            let suffix_len = rng.below(80) as usize;
            let subject_text: Vec<u8> = (0..prefix_len + core_len + suffix_len)
                .map(|_| {
                    if rng.below(1000) < subject_ambiguity {
                        AMBIGUOUS[rng.below(AMBIGUOUS.len() as u64) as usize]
                    } else {
                        BASES[rng.below(BASES.len() as u64) as usize]
                    }
                })
                .collect();
            let len_full = subject_text.len();
            let packed = encode_subject_ncbi2na_packed(&subject_text);
            let codes = resolve_ncbi4na_to_ncbi2na(
                &subject_text
                    .iter()
                    .map(|&letter| iupacna_to_ncbi4na(letter).unwrap_or(15))
                    .collect::<Vec<u8>>(),
            );
            let unpacked_direct = x_unpack_ncbi2na(&packed, len_full);
            assert_eq!(unpacked_direct, codes, "case {case}: unpacked copy");

            // query (blastna): flank, the subject core with changes, flank
            let mut query: Vec<u8> = Vec::new();
            for _ in 0..rng.below(40) {
                query.push(rng.below(4) as u8);
            }
            let mut anchors: Vec<(usize, usize)> = Vec::new();
            for (s, &base) in codes.iter().enumerate().skip(prefix_len).take(core_len) {
                let roll = rng.below(1000);
                if roll < indel {
                    // extra base in the subject
                } else if roll < 2 * indel {
                    // extra base in the query
                    query.push(rng.below(4) as u8);
                    anchors.push((query.len(), s));
                    query.push(base);
                } else if roll < 2 * indel + substitution {
                    query.push((base + 1 + rng.below(3) as u8) & 3);
                } else if roll < 2 * indel + substitution + query_ambiguity {
                    // an ambiguity code in the query (never matches)
                    query.push(4 + rng.below(12) as u8);
                } else {
                    anchors.push((query.len(), s));
                    query.push(base);
                }
            }
            for _ in 0..rng.below(40) {
                query.push(rng.below(4) as u8);
            }

            // the chunk: starts at a multiple of 4 after the subject start (sometimes inside
            // the core), and is sometimes shorter than the rest of the subject
            let start_limit = if rng.below(4) == 0 {
                prefix_len + core_len / 2
            } else {
                prefix_len
            };
            let chunk_offset = 4 * (1 + rng.below((start_limit / 4) as u64) as usize);
            let rest = len_full - chunk_offset;
            let chunk_len = if rng.below(3) == 0 {
                rest
            } else {
                1 + rng.below(rest as u64) as usize
            };
            let packed_chunk = &packed[chunk_offset / COMPRESSION_RATIO..];
            let x_subject = XPrelimSubject::new(true, &packed, len_full);
            let unpacked_chunk = x_subject.chunk(chunk_offset).expect("switch on");
            assert_eq!(unpacked_chunk, &codes[chunk_offset..], "case {case}: chunk");

            // one seed per phase s_off % 4 (when the chunk has one)
            let scoring = match case % 5 {
                0 => (2, -3, 5, 2),
                1 => (1, -3, 0, 0),
                _ => (1, -2, 0, 0),
            };
            let x_drop = 5 + rng.below(60) as i32;
            for phase in 0..4 {
                let candidates: Vec<(usize, usize)> = anchors
                    .iter()
                    .filter(|&&(_, s)| {
                        s >= chunk_offset && s < chunk_offset + chunk_len && s % 4 == phase
                    })
                    .map(|&(q, s)| (q, s - chunk_offset))
                    .collect();
                if candidates.is_empty() {
                    continue;
                }
                let (q_off, s_off) = candidates[rng.below(candidates.len() as u64) as usize];
                let run = |subject: Option<&[u8]>, scratch: &mut GreedyAlignScratch| {
                    x_greedy_gapped_alignment_score_only(
                        subject,
                        &query,
                        packed_chunk,
                        chunk_len,
                        q_off,
                        s_off,
                        scoring.0,
                        scoring.1,
                        scoring.2,
                        scoring.3,
                        x_drop,
                        scratch,
                    )
                };
                // the packed call as ported
                let reference = run(None, &mut GreedyAlignScratch::new());
                let fresh = run(Some(unpacked_chunk), &mut GreedyAlignScratch::new());
                assert_eq!(fresh, reference, "case {case}: phase {phase}, new scratch");
                let reused = run(Some(unpacked_chunk), &mut used);
                assert_eq!(
                    reused, reference,
                    "case {case}: phase {phase}, used scratch"
                );
                calls += 1;
                phases[phase] += 1;
                if let Some(result) = reference {
                    aligned += 1;
                    if result.1 - result.0 >= 64 {
                        long_runs += 1;
                    }
                }
            }
        }
        // the generator must produce real alignments at every phase
        assert!(
            aligned * 2 > calls,
            "only {aligned} alignments in {calls} calls"
        );
        assert!(
            phases.iter().all(|&n| n * 8 > cases),
            "seeds per phase {phases:?} in {cases} cases"
        );
        assert!(long_runs * 4 > calls, "only {long_runs} long alignments");
    }
}
