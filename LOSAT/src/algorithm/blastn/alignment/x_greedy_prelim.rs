//! No NCBI counterpart: the preliminary (score-only) greedy gapped alignment of megablast run on
//! the packed ncbi2na subject widened to one base per byte, with `compressed_subject = false`, so
//! that the uncompressed scans of `BLAST_GreedyAlign` serve the preliminary stage too. A
//! score-only call on those bytes runs `x_greedy_fused`, and `x_greedy_row` only where the fused
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
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:516-520
//! ```c
//!           if (fence_hit) {
//!               hsp_list->hsp_array[index] = Blast_HSPFree(hsp_list->hsp_array[index]);
//!               break;
//!           }
//! ```
//! A call widens only the bases around its ungapped hit (`X_WINDOW_MARGIN` on each side, into
//! the worker's `GreedyAlignScratch::x_window`), with NCBI's fence sentry before and after them
//! where the window does not reach an end of the chunk, as NCBI's partial fetch places it
//! (seqdbvol.cpp:1565-1595, `x_nucleotide_blastna_ranges`); the next call of the same chunk on
//! that worker reads the same window when its hit and margins lie inside it (a chunk shorter than
//! the margins is widened once per worker). The scans test the sentry where a
//! slide stops (greedy_align.c), so a call that does not meet it read only bases of the window,
//! the same bytes as the whole copy holds, at the same offsets less the window's start, and its
//! result is the whole copy's. No state of the scans starts past a sentry without one meeting
//! it: a slide that starts at the sentry reads it, and an affine deletion run, which passes
//! subject bases without reading them, reaches the sentry's offset or beyond only in states the
//! X-drop test drops (each further base lowers the run's score by the extension cost while the
//! best score does not fall); with no sentry on a side, that side's length is the whole copy's. A call that meets it returns no result (as NCBI's traceback drops
//! the HSP) and runs again on the one-byte copy of the whole subject, made on the first such
//! call and shared by every chunk and worker of that subject, from the `max_dist` it started
//! with (the fenced call may have doubled it). `max_dist` is sized by the chunk's length in both
//! calls (`greedy_gapped_alignment_internal`).

use super::greedy::{greedy_gapped_alignment_internal, GreedyAlignScratch, FENCE_SENTRY};
use super::x_fence;
use crate::core::blast_encoding::COMPRESSION_RATIO;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: the bases on each side of its ungapped hit that a preliminary call's
/// window holds. A size, not a threshold: a call that reads past it meets a fence and runs again
/// on the whole subject (module doc), so the result does not depend on it.
const X_WINDOW_MARGIN: usize = 1024;

/// No NCBI counterpart: the last number given to an `XPrelimSubject` that made a window, so that
/// a worker's window names the subject it holds.
static X_LAST_SUBJECT: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: a worker's window (`XPrelimSubject::score_only`): bases `from..to` of the
/// chunk that starts at base `chunk_offset` of subject number `subject` (0: none), one per byte,
/// at `bytes[lead..lead + to - from]`, between fence sentries; `len` bytes in all. `bytes` only
/// grows; bytes past `len` are left from earlier windows. The subject and the chunk's offset fix
/// the chunk's length (both callers pass the chunk's `length`), so they fix where the sentries
/// stand.
#[derive(Default)]
pub(crate) struct XPrelimWindow {
    bytes: Vec<u8>,
    subject: u64,
    chunk_offset: usize,
    from: usize,
    to: usize,
    lead: usize,
    len: usize,
}

// No NCBI counterpart (memory placement only): a window past the steady size goes back to none
// (`scratch_pool::Settle`); the next call fills a window as a new scratch's first call does.
impl crate::utils::scratch_pool::Settle for XPrelimWindow {
    fn settle(&mut self) {
        if crate::utils::scratch_pool::passes(crate::utils::scratch_pool::vec_bytes(&self.bytes)) {
            *self = Self::default();
        }
    }
}

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
/// One subject's packed sequence for the preliminary greedy calls: windows of it (`score_only`)
/// and the one-byte copy of all of it, made on the first call that meets a fence and shared by
/// every chunk and worker of that subject. Off unless this is the preliminary stage of a search
/// that uses the greedy aligner. No NCBI counterpart for the windows and the copy: NCBI reads the
/// packed bytes.
pub(crate) struct XPrelimSubject<'a> {
    packed_full: &'a [u8],
    len_full: usize,
    enabled: bool,
    unpacked: OnceLock<Vec<u8>>,
    number: OnceLock<u64>,
}

impl<'a> XPrelimSubject<'a> {
    /// `greedy_prelim`: the preliminary stage of a search that uses the greedy aligner.
    pub(crate) fn new(greedy_prelim: bool, packed_full: &'a [u8], len_full: usize) -> Self {
        Self {
            packed_full,
            len_full,
            enabled: greedy_prelim,
            unpacked: OnceLock::new(),
            number: OnceLock::new(),
        }
    }

    /// The copy from base `chunk_offset` (a multiple of 4) to the end of the subject: the bases
    /// of `packed_full[chunk_offset / 4..]`. `None` when the copy is off.
    fn chunk(&self, chunk_offset: usize) -> Option<&[u8]> {
        if !self.enabled {
            return None;
        }
        debug_assert_eq!(chunk_offset % COMPRESSION_RATIO, 0);
        let unpacked = self
            .unpacked
            .get_or_init(|| x_unpack_ncbi2na(self.packed_full, self.len_full));
        unpacked.get(chunk_offset..)
    }

    /// Writes bases `from..to` of the chunk that starts at base `chunk_offset` (the base
    /// `chunk_offset + from` a multiple of 4, `to <= subject_len`, the chunk's length) into
    /// `window`, one per byte, after `FENCE_SENTRY` when `from > 0` and followed by it when
    /// `to < subject_len`.
    fn x_fill_window(
        &self,
        window: &mut XPrelimWindow,
        chunk_offset: usize,
        from: usize,
        to: usize,
        subject_len: usize,
    ) {
        debug_assert_eq!((chunk_offset + from) % COMPRESSION_RATIO, 0);
        let lead = usize::from(from > 0);
        let first = (chunk_offset + from) / COMPRESSION_RATIO;
        let bytes = (to - from).div_ceil(COMPRESSION_RATIO);
        let needed = lead + bytes * COMPRESSION_RATIO + 1;
        if window.bytes.len() < needed {
            window.bytes.resize(needed, 0);
        }
        if lead == 1 {
            window.bytes[0] = FENCE_SENTRY;
        }
        for (bases, &byte) in window.bytes[lead..lead + bytes * COMPRESSION_RATIO]
            .chunks_exact_mut(COMPRESSION_RATIO)
            .zip(&self.packed_full[first..first + bytes])
        {
            bases.copy_from_slice(&X_UNPACKED[byte as usize]);
        }
        let mut len = lead + (to - from);
        if to < subject_len {
            window.bytes[len] = FENCE_SENTRY;
            len += 1;
        }
        window.subject = *self
            .number
            .get_or_init(|| X_LAST_SUBJECT.fetch_add(1, Ordering::Relaxed) + 1);
        window.chunk_offset = chunk_offset;
        window.from = from;
        window.to = to;
        window.lead = lead;
        window.len = len;
    }

    /// The preliminary greedy alignment (`x_greedy_gapped_alignment_score_only`) seeded at
    /// `(q_off, s_off)` of the chunk that starts at base `chunk_offset` (a multiple of 4) and has
    /// `subject_len` bases, for the ungapped hit at `hit` (subject offsets of the chunk): first on
    /// a window that holds the hit and `X_WINDOW_MARGIN` bases on each side, again on the whole
    /// copy when that call meets a fence (module doc). Subject offsets of the result are offsets
    /// of the chunk.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn score_only(
        &self,
        chunk_offset: usize,
        subject_len: usize,
        hit: (usize, usize),
        query: &[u8],
        q_off: usize,
        s_off: usize,
        reward: i32,
        penalty: i32,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        scratch: &mut GreedyAlignScratch,
    ) -> Option<(usize, usize, usize, usize, i32, usize, usize)> {
        assert!(
            self.enabled,
            "a score-only greedy call runs only in the preliminary stage of a greedy search"
        );
        if s_off >= subject_len {
            // `greedy_gapped_alignment_internal` gives no result (no subject base to the right).
            return None;
        }
        let from = hit.0.min(s_off).saturating_sub(X_WINDOW_MARGIN) & !(COMPRESSION_RATIO - 1);
        let to = (hit.1.max(s_off + 1) + X_WINDOW_MARGIN).min(subject_len);
        let mut window = std::mem::take(&mut scratch.x_window);
        let held = self.number.get().is_some_and(|&n| n == window.subject)
            && window.chunk_offset == chunk_offset
            && window.from <= from
            && to <= window.to;
        if !held {
            self.x_fill_window(&mut window, chunk_offset, from, to, subject_len);
        }
        let max_dist = scratch.max_dist;
        let fences = x_fence::x_fence_count();
        let core = greedy_gapped_alignment_internal(
            query,
            &window.bytes[..window.len],
            window.len,
            subject_len,
            q_off,
            s_off - window.from + window.lead,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
            false,
        );
        let shift = window.from - window.lead;
        scratch.x_window = window;
        if x_fence::x_fence_count() == fences {
            return x_score_only_result(core, shift);
        }
        scratch.max_dist = max_dist;
        x_greedy_gapped_alignment_score_only(
            self.chunk(chunk_offset)
                .expect("the copy is on wherever a score-only greedy call runs"),
            query,
            subject_len,
            q_off,
            s_off,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
        )
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
/// same result as the packed call). LOSAT has no packed form of this call. `score_only` runs it
/// on a window first.
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
    x_score_only_result(core, 0)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2539-2569
/// (s_BlastGreedyGapAlignStructFill: query_start/stop, subject_start/stop, score)
/// The result of a score-only call whose subject bytes start `shift` bases into the chunk (a
/// window's start less its sentry, 0 for the whole copy), in offsets of the chunk. A window call
/// that did not meet a fence starts at or after the window's first base, so the shift never
/// turns a negative start into a valid one.
fn x_score_only_result(
    core: Option<super::greedy::GreedyGappedCore>,
    shift: usize,
) -> Option<(usize, usize, usize, usize, i32, usize, usize)> {
    let core = core?;
    if core.q_start < 0 || core.s_start < 0 {
        return None;
    }
    Some((
        core.q_start as usize,
        core.q_end as usize,
        core.s_start as usize + shift,
        core.s_end as usize + shift,
        core.score,
        core.q_seed_start as usize,
        core.s_seed_start as usize + shift,
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

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:516-520
    /// ```c
    ///           if (fence_hit) {
    ///               hsp_list->hsp_array[index] = Blast_HSPFree(hsp_list->hsp_array[index]);
    ///               break;
    ///           }
    /// ```
    /// A sequence of score-only calls through `XPrelimSubject::score_only` (windows kept across
    /// calls, the whole copy after a fence) gives each call's result on the whole copy and leaves
    /// `max_dist` where the whole-copy calls leave it: hits of a query made of mutated pieces of
    /// the subject (alignments within the margins and past them), chunks at offset 0 and inside
    /// the subject, hits at the ends of chunks, linear and affine scores.
    #[test]
    fn score_only_on_windows_gives_the_whole_copy_result() {
        const SUBJECT_LEN: usize = 30_000;
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let subject: Vec<u8> = (0..SUBJECT_LEN).map(|_| (next() % 4) as u8).collect();
        let text: Vec<u8> = subject.iter().map(|&code| b"ACGT"[code as usize]).collect();
        let packed = encode_subject_ncbi2na_packed(&text);
        // Pieces of the subject, 1-5 % substituted, some with an insertion or a deletion; 12
        // seeds spread over each piece, at the piece's offsets before its indel.
        let mut query = Vec::new();
        let mut seeds = Vec::new();
        for (piece, &(start, len)) in [
            (0usize, 300usize),
            (2_900, 6_000),
            (11_500, 900),
            (11_990, 40),
            (15_000, 2_500),
            (19_900, 200),
            (24_000, 5_999),
            (29_700, 300),
        ]
        .iter()
        .enumerate()
        {
            let q_start = query.len();
            let substitution = 100 / (piece as u64 % 5 + 1);
            let indel_at = if piece % 3 == 1 { len / 2 } else { len };
            for i in 0..len {
                if i == indel_at && piece % 2 == 0 {
                    continue;
                }
                if i == indel_at {
                    query.push((next() % 4) as u8);
                }
                let base = subject[start + i];
                query.push(if next() % substitution == 0 {
                    (base + 1 + (next() % 3) as u8) % 4
                } else {
                    base
                });
            }
            query.extend((0..50).map(|_| (next() % 4) as u8));
            let span = indel_at.min(len);
            seeds.extend((0..12).map(|k| (q_start + k * span / 12, start + k * span / 12)));
        }
        let whole = x_unpack_ncbi2na(&packed, SUBJECT_LEN);
        let (mut calls, mut fenced) = (0usize, 0usize);
        for (chunk_offset, subject_len) in [
            (0usize, SUBJECT_LEN),
            (0, 12_000),
            (8_000, 12_000),
            (11_996, 18_004),
        ] {
            for (reward, penalty, gap_open, gap_extend, x_drop) in
                [(1, -2, 0, 0, 20), (1, -3, 0, 0, 40), (2, -3, 5, 2, 30)]
            {
                let subject = XPrelimSubject::new(true, &packed, SUBJECT_LEN);
                let (mut windows, mut copy) =
                    (GreedyAlignScratch::new(), GreedyAlignScratch::new());
                for &(q_off, s_seed) in &seeds {
                    if s_seed < chunk_offset || s_seed >= chunk_offset + subject_len {
                        continue;
                    }
                    let s_off = s_seed - chunk_offset;
                    let half = (next() % 400) as usize;
                    let hit = (
                        s_off.saturating_sub(half),
                        (s_off + half).min(subject_len - 1),
                    );
                    let fences = x_fence::x_fence_count();
                    let got = subject.score_only(
                        chunk_offset,
                        subject_len,
                        hit,
                        &query,
                        q_off,
                        s_off,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        x_drop,
                        &mut windows,
                    );
                    fenced += usize::from(x_fence::x_fence_count() != fences);
                    calls += 1;
                    let want = x_greedy_gapped_alignment_score_only(
                        &whole[chunk_offset..],
                        &query,
                        subject_len,
                        q_off,
                        s_off,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        x_drop,
                        &mut copy,
                    );
                    assert_eq!(
                        got, want,
                        "chunk {chunk_offset}+{subject_len}, scores {reward}/{penalty} \
                         {gap_open}/{gap_extend}, seed ({q_off}, {s_off})"
                    );
                    assert_eq!(windows.max_dist, copy.max_dist);
                }
            }
        }
        assert!(
            fenced * 20 > calls && fenced * 2 < calls,
            "calls {calls}, fenced {fenced}"
        );
    }
}
