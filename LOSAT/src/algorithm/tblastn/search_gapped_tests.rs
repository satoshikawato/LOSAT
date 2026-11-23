use super::*;
use crate::algorithm::blastp::gapalign::blast_gapped_alignment_with_traceback_with_scratch;
use crate::algorithm::tblastn::stage_d_linking::{
    link_preliminary_hsps, reap_by_evalue, LinkedHspList,
};
use crate::algorithm::tblastn::stage_d_stats::{
    local_parameters_for_call, LocalParameterCall, LocalParameterOptions,
};
use crate::config::ProteinScoringSpec;
use crate::stats::spouge::lookup_protein_gumbel_params;
use crate::stats::tables::{lookup_protein_params_gapped, lookup_protein_params_ungapped};
use std::fs;

/// A frame of `letters` (NCBISTDAA) with its two NULLBs.
fn one_hit_frame(frame: i8, letters: &[u8]) -> QueryFrame {
    let mut aa_seq = vec![0];
    aa_seq.extend_from_slice(letters);
    aa_seq.push(0);
    QueryFrame {
        frame,
        aa_seq,
        aa_seq_nomask: None,
        aa_len: letters.len(),
        orig_len: 3 * letters.len(),
        seg_masks: Vec::new(),
    }
}

/// The PSSM of line RG item X2d's probe (without its motif rows; the fixtures' `x2d_d50.asn`): 30
/// rows of -6 ('-' BLAST_SCORE_MIN), row 9 -61 for every letter, W +60 at row 10, V -60 at row
/// 11, K +1 at row 12, S +4 at row 14 and '-' `dash` at rows 13 to 15.
fn x2d_pssm(dash: i32) -> crate::algorithm::pssm_search::rows::QueryPssm {
    const W: usize = 20;
    const V: usize = 19;
    const K: usize = 10;
    const S: usize = 17;
    let rows = (0..30)
        .map(|r| {
            let mut row = [-6; crate::core::blast_util::BLASTAA_SIZE];
            if r == 9 {
                row = [-61; crate::core::blast_util::BLASTAA_SIZE];
            }
            row[0] = -32768;
            match r {
                10 => row[W] = 60,
                11 => row[V] = -60,
                12 => row[K] = 1,
                14 => row[S] = 4,
                _ => {}
            }
            if (13..=15).contains(&r) {
                row[0] = dash;
            }
            row
        })
        .collect();
    crate::algorithm::pssm_search::rows::QueryPssm { rows }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4250-4258,4295,
// 4342-4345 (cited at `x_score_only_past_chunk`) and :3405-3439
// The word W V K at letters 14..17 of a frame of 19 letters (its last two N S), with the
// probe's PSSM: a one-hit HSP of length -1 at (10, 14) whose window of 11 rows reads the
// frame's last five letters, its NULLB (row 15, '-' 50) and five letters of the next frame
// (A, -6): 60 - 60 + 1 - 6 + 4 + 50 - 30 = 19 > 0, so the gapped start is (15, 19), on the
// NULLB. NCBI's left extension reads back from the NULLB: 50, S/S, three letters of the query
// and three of the subject in gaps (-14 each), W/W: 86 at (10, 14), the score of NCBI 2.17.0's
// alignment there (37.6 bits, raw 86, in the search of `db.x2d.tblastn.nullb.db.6` at `-evalue
// 10`, which LOSAT refuses at the traceback); the right extension does not run, so the HSP
// ends at the start.
#[test]
fn a_one_hit_gapped_start_on_its_frames_nullb_reads_it_with_the_pssms_dash_column() {
    const A: u8 = 1;
    let pssm = x2d_pssm(50);
    let mut letters = vec![A; 14];
    letters.extend_from_slice(&[20, 19, 10, 13, 17]);
    let frame = one_hit_frame(1, &letters);
    let query_frame = one_hit_frame(0, &[A; 30]);
    let hit = InitHsp {
        frame: 1,
        chunk_offset: 0,
        q_seed: 10,
        s_seed: 14,
        q_start: 10,
        s_start: 14,
        length: -1,
        score: 60,
    };
    let query_frames = [query_frame];
    // The six frames' buffer of a subject whose frames hold `first` letters, then one fewer.
    let buffers = |first: usize, after_frame: Vec<u8>| {
        let mut frame_offsets = [0; 7];
        for k in 1..7 {
            frame_offsets[k] = frame_offsets[k - 1] + first - usize::from(k > 1) + 1;
        }
        XOneHitBuffers {
            frame_offsets,
            frame_index: 0,
            after_frame,
            frame_end: std::cell::OnceCell::new(),
        }
    };
    let gapped = |pssm: &crate::algorithm::pssm_search::rows::QueryPssm,
                  frame: &QueryFrame,
                  buffers: &XOneHitBuffers| {
        let chunk = TranslatedChunk::new(0, frame.aa_len);
        score_gapped_chunk::<_, true>(
            &[&[A; 30]],
            &query_frames,
            &[0],
            30,
            frame,
            chunk,
            &frame.aa_seq[1..=frame.aa_len],
            &frame.aa_seq[1..],
            &[hit],
            11,
            1,
            38,
            &[1],
            ScoringMatrix::Blosum62,
            &mut (
                GapAlignScratch::new(),
                crate::algorithm::pssm_search::gapalign::PssmGapScratch::new(),
            ),
            super::super::score::PssmScore(pssm),
            Some(buffers),
        )
        .unwrap()
    };
    let hsp = |score| {
        vec![(
            0,
            GappedHsp {
                frame: 1,
                score,
                q_start: 10,
                q_end: 15,
                q_gapped_start: 15,
                s_start: 14,
                s_end: 19,
                s_gapped_start: 19,
            },
        )]
    };
    assert_eq!(gapped(&pssm, &frame, &buffers(19, vec![A; 10])), hsp(86));
    // Without the frame's last letter S the word ends one letter before the NULLB (row 14) and
    // the start (15, 19) lies on the next frame's first letter: the window is 60 - 60 + 1 - 6
    // + 50 - 36 = 9, and the left extension reads that letter (A -6, or with row 15 scoring A
    // +9 instead: 9), the NULLB (50) and the frame's letters back to W/W over the gaps (-14
    // each): 76 (91) at (10, 14), up to the start; a W there scores -6 in both.
    let frame = one_hit_frame(1, &letters[..18]);
    let mut row15 = x2d_pssm(50);
    row15.rows[15][usize::from(A)] = 9;
    for (pssm, first, score) in [(&pssm, A, 76), (&row15, A, 91), (&row15, 20, 76)] {
        let mut after = vec![A; 10];
        after[0] = first;
        assert_eq!(gapped(pssm, &frame, &buffers(18, after)), hsp(score));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1058,1063,1098-1101 and
// blast_gapalign.c:4251-4258 (cited at `x_gapped_start_past_translation`)
// The word W V K at letters 14..17 of frame -3, the last of the translation buffer, with the
// probe's PSSM whose '-' column scores 50 at rows 13 to 20: the one-hit HSP of length -1 at
// (10, 14) has a window of 11 rows that reads the frame's last letters, the NULLB after frame -3
// and NULLBs past NCBI's buffer (exception 3), which sum above 0, so its gapped start is (15,
// 19). With 19 letters that start is the NULLB after frame -3, the buffer's last byte: the
// search goes on. With 18 letters it lies past that NULLB, where NCBI's left extension reads a
// byte the buffer never sets: refused (rule 3).
#[test]
fn a_one_hit_gapped_start_past_the_translation_buffer_is_refused() {
    const A: u8 = 1;
    let mut pssm = x2d_pssm(50);
    for row in &mut pssm.rows[16..=20] {
        row[0] = 50;
    }
    let mut letters = vec![A; 14];
    letters.extend_from_slice(&[20, 19, 10, 13, 17]);
    let query_frames = [one_hit_frame(0, &[A; 30])];
    let hit = InitHsp {
        frame: -3,
        chunk_offset: 0,
        q_seed: 10,
        s_seed: 14,
        q_start: 10,
        s_start: 14,
        length: -1,
        score: 60,
    };
    let gapped = |aa_len: usize| {
        let frame = one_hit_frame(-3, &letters[..aa_len]);
        // Five frames of 20 letters before it, each after its NULLB.
        let mut frame_offsets = [0; 7];
        for k in 1..6 {
            frame_offsets[k] = frame_offsets[k - 1] + 21;
        }
        frame_offsets[6] = frame_offsets[5] + aa_len + 1;
        let buffers = XOneHitBuffers {
            frame_offsets,
            frame_index: 5,
            after_frame: Vec::new(),
            frame_end: std::cell::OnceCell::new(),
        };
        let chunk = TranslatedChunk::new(0, aa_len);
        score_gapped_chunk::<_, true>(
            &[&[A; 30]],
            &query_frames,
            &[0],
            30,
            &frame,
            chunk,
            &frame.aa_seq[1..=aa_len],
            &frame.aa_seq[1..],
            &[hit],
            11,
            1,
            38,
            &[1],
            ScoringMatrix::Blosum62,
            &mut (
                GapAlignScratch::new(),
                crate::algorithm::pssm_search::gapalign::PssmGapScratch::new(),
            ),
            super::super::score::PssmScore(&pssm),
            Some(&buffers),
        )
    };
    let on_nullb = gapped(19).unwrap();
    assert!(
        !on_nullb.is_empty() && on_nullb.iter().all(|(_, hsp)| hsp.s_gapped_start == 19),
        "{on_nullb:?}"
    );
    let message = gapped(18).unwrap_err().to_string();
    assert!(
        message.contains("lies past the end of the subject's translation"),
        "{message}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4609 and :569-575
// (cited at `x_traceback_start_past_frame`)
/// The traceback's check of a gapped start against its frame's NULLB (`aa_len`, here 6 for
/// frame +1 of 19 bases): before it and on it the traceback runs, after it NCBI reads a byte its
/// translation does not set (rule 3), in frame coordinates whatever the window (`subject_base`
/// -1 for a window from the first letter with its left fence, 1 for one from letter 2).
#[test]
fn the_traceback_refuses_a_gapped_start_after_its_frames_nullb() {
    for (s_start, base, refused) in [(5, 0, false), (6, 0, false), (7, -1, false), (7, 0, true)] {
        let result = x_traceback_start_in_frame(19, 1, s_start, base);
        assert_eq!(result.is_err(), refused, "{s_start} {base}");
    }
    let message = x_traceback_start_in_frame(19, 1, 6, 1)
        .unwrap_err()
        .to_string();
    assert!(
        message.contains("starts after the NULLB after its frame"),
        "{message}"
    );
    // A subject of fewer than 3 bases has no frame +1: an internal error.
    assert!(x_traceback_start_in_frame(2, 1, 0, 0).is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4171-4189,4257-4258,
// 4281,4295,4342-4345 (cited at `x_score_only_past_chunk`)
// A one-hit HSP of length -1 at letter `chunk.length - 3` of a chunk that is not its frame's
// last, whose window of 11 letters (all of them past the chunk's end but 3) sums above 0: its
// gapped start lies 2 letters past the chunk (`q_start + 5`). NCBI's left extension reads back
// from it over the frame's letters (`AdjustSubjectRange` leaves a chunk shorter than
// MAX_SUBJECT_OFFSET as it is), and the right one does not run (`s_length` is past
// `subject->length`): the HSP is the 11 letters W/W that end at the gapped start (11 x 11),
// although the letters W/W go on to the right.
#[test]
fn a_one_hit_gapped_start_past_its_chunk_extends_only_to_the_left() {
    const W: u8 = 20;
    const A: u8 = 1;
    let query = [W; 20];
    let mut letters = vec![A; 60];
    letters[22..50].fill(W);
    let frame = one_hit_frame(1, &letters);
    let query_frame = one_hit_frame(0, &query);
    let chunk = TranslatedChunk::new(0, 30);
    let hit = InitHsp {
        frame: 1,
        chunk_offset: 0,
        q_seed: 5,
        s_seed: 27,
        q_start: 5,
        s_start: 27,
        length: -1,
        score: 4,
    };
    let query_frames = [query_frame];
    let buffers = XOneHitBuffers {
        frame_offsets: [0, 61, 121, 181, 242, 302, 362],
        frame_index: 0,
        after_frame: Vec::new(),
        frame_end: std::cell::OnceCell::new(),
    };
    let gapped = score_gapped_chunk::<_, true>(
        &[&query],
        &query_frames,
        &[0],
        20,
        &frame,
        chunk,
        &frame.aa_seq[1..31],
        &frame.aa_seq[1..],
        &[hit],
        11,
        1,
        38,
        &[1],
        ScoringMatrix::Blosum62,
        &mut GapAlignScratch::new(),
        super::super::score::MatrixScore,
        Some(&buffers),
    )
    .unwrap();
    assert_eq!(
        gapped,
        vec![(
            0,
            GappedHsp {
                frame: 1,
                score: 121,
                q_start: 0,
                q_end: 10,
                q_gapped_start: 10,
                s_start: 22,
                s_end: 32,
                s_gapped_start: 32,
            }
        )]
    );
    // The same HSP without the buffers of the one-hit word finder is an internal error.
    let error = score_gapped_chunk::<_, true>(
        &[&query],
        &[one_hit_frame(0, &query)],
        &[0],
        20,
        &frame,
        chunk,
        &frame.aa_seq[1..31],
        &frame.aa_seq[1..],
        &[hit],
        11,
        1,
        38,
        &[1],
        ScoringMatrix::Blosum62,
        &mut GapAlignScratch::new(),
        super::super::score::MatrixScore,
        None,
    )
    .unwrap_err();
    assert!(error.to_string().contains("internal error"), "{error}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1058,1098-1101 (cited at
// `XOneHitBuffers`)
// The subject's window reads NCBI's bytes: past a frame its NULLB and the next frame's letters;
// past the NULLB after frame -3 every letter reads as the sentinel NULLB (PD-LOSAT-NCBI-DEFECTS
// exception 3). The query's window is the context's letters and NULLBs past them
// (`gapped_start_window`).
#[test]
fn one_hit_windows_read_ncbi_buffers_and_the_sentinel_past_them() {
    // Frames of 5, 4, 4, 5, 4, 4 letters (a subject of 15 bases): NULLBs at 0, 6, 11, 16, 22,
    // 27 and 32.
    let frames: Vec<QueryFrame> = [(1, 5), (2, 4), (3, 4), (-1, 5), (-2, 4), (-3, 4)]
        .into_iter()
        .enumerate()
        .map(|(index, (frame, length))| one_hit_frame(frame, &vec![10 + index as u8; length]))
        .collect();
    let whole = |frame_index: usize| TranslatedChunk::new(0, frames[frame_index].aa_len);
    let buffers = |frame_index: usize| {
        XOneHitBuffers::new(
            &frames,
            frame_index,
            whole(frame_index),
            super::super::x_subject_store::XPrelimSubject::Nt(SubjectNt::Letters(&[])),
            1,
        )
        .unwrap()
    };
    let first = buffers(0);
    assert_eq!(first.frame_offsets, [0, 6, 11, 16, 22, 27, 32]);
    assert_eq!(
        first.after_frame,
        vec![11, 11, 11, 11, 0, 12, 12, 12, 12, 0]
    );
    let window = |buffers: &XOneHitBuffers, frame_index: usize, s: usize| {
        let frame = &frames[frame_index];
        buffers.subject_window(frame, whole(frame_index), &frame.aa_seq[1..], s)
    };
    // A context's window: its letters, then NULLBs (the first one is the context's own).
    assert_eq!(
        gapped_start_window(&[1, 2, 3, 4], 0),
        [1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0]
    );
    assert_eq!(
        gapped_start_window(&[5, 6, 7, 8, 9, 10, 11, 12], 7),
        [12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    );
    // Frame 1's letter 0: the window holds frame 1's NULLB and 5 bytes of frame 2.
    let s = window(&first, 0, 0).unwrap();
    assert_eq!(s, [10, 10, 10, 10, 10, 0, 11, 11, 11, 11, 0]);
    // Frame -1 (NULLB at 16): its letters 0..10 are bytes 17..27, the last one frame -2's
    // NULLB; from its letter 1 the window ends on byte 28, frame -3's first letter.
    let fourth = buffers(3);
    let s = window(&fourth, 3, 0).unwrap();
    assert_eq!(s, [13, 13, 13, 13, 13, 0, 14, 14, 14, 14, 0]);
    let s = window(&fourth, 3, 1).unwrap();
    assert_eq!(s, [13, 13, 13, 13, 0, 14, 14, 14, 14, 0, 15]);
    // Frame -2 (NULLB at 22): from its letter 0 (byte 23) the window ends on byte 33, past
    // the NULLB after frame -3 (byte 32), which reads as a NULLB; from frame -3's letters the
    // window holds its NULLB and then NULLBs.
    let fifth = buffers(4);
    assert_eq!(fifth.after_frame, vec![15, 15, 15, 15, 0]);
    let s = window(&fifth, 4, 0).unwrap();
    assert_eq!(s, [14, 14, 14, 14, 0, 15, 15, 15, 15, 0, 0]);
    let sixth = buffers(5);
    assert!(sixth.after_frame.is_empty());
    let s = window(&sixth, 5, 0).unwrap();
    assert_eq!(s, [15, 15, 15, 15, 0, 0, 0, 0, 0, 0, 0]);
    let s = window(&sixth, 5, 3).unwrap();
    assert_eq!(s, [15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    assert!(sixth.past_translation(5) && !sixth.past_translation(4));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250 and
// c++/src/algo/blast/core/blast_util.c:1098-1101 (cited at `XOneHitBuffers`)
// A chunk that is not its frame's last can end fewer than HSP_MAX_WINDOW letters before the
// frame does (the frame goes on 1 letter past it): a window from its end reads the frame's last
// letters, its NULLB and the next frame's letters, as from the frame's last chunk. A next frame
// longer than the window is cut without its NULLB (NCBI's byte there is a letter).
#[test]
fn one_hit_windows_of_a_chunk_near_its_frame_end_read_the_next_frame() {
    let frames: Vec<QueryFrame> = [(1, 31), (2, 31), (3, 30), (-1, 31), (-2, 31), (-3, 30)]
        .into_iter()
        .enumerate()
        .map(|(index, (frame, length))| one_hit_frame(frame, &vec![10 + index as u8; length]))
        .collect();
    let chunk = TranslatedChunk::new(0, 30);
    let buffers = XOneHitBuffers::new(
        &frames,
        0,
        chunk,
        super::super::x_subject_store::XPrelimSubject::Nt(SubjectNt::Letters(&[])),
        1,
    )
    .unwrap();
    assert_eq!(buffers.after_frame, vec![11; 10]);
    let frame = &frames[0];
    let s = buffers
        .subject_window(frame, chunk, &frame.aa_seq[1..], 26)
        .unwrap();
    assert_eq!(s, [10, 10, 10, 10, 10, 0, 11, 11, 11, 11, 11]);
    // A chunk that ends HSP_MAX_WINDOW letters or more before its frame's end holds nothing
    // after the frame.
    let early = TranslatedChunk::new(0, 21);
    let buffers = XOneHitBuffers::new(
        &frames,
        0,
        early,
        super::super::x_subject_store::XPrelimSubject::Nt(SubjectNt::Letters(&[])),
        1,
    )
    .unwrap();
    assert!(buffers.after_frame.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1063-1064,1080,1098-1101
// (cited at `XOneHitBuffers`)
// The bytes after a frame of a plan without residues (LOSAT_X_TBNUNITXLATE, the large subjects'
// plans), which the unit translates itself, are those of the plan with residues: the next
// frames' letters, each with its NULLB, as NCBI's translation buffer holds them; for every
// frame, for its last chunk and the chunks that end fewer than HSP_MAX_WINDOW letters before it.
#[test]
fn one_hit_bytes_after_a_unit_translated_frame_are_the_translation_buffer() {
    use super::super::x_subject_store::XPrelimSubject;
    use super::super::x_subject_translate::{x_packed_lean_plan, x_packed_subject_plan};
    use crate::sequence::packed_subject::PackedNucleotide;
    for length in [14usize, 20, 31, 47, 302, 3001] {
        let letters: Vec<u8> = (0..length)
            .map(|i| b"ACGT"[(i * 7 + i / 3 + i / 11) % 4])
            .collect();
        let store = PackedNucleotide::from_subject_letters(&letters).unwrap();
        let full = x_packed_subject_plan(&store, 1).unwrap();
        let lean = x_packed_lean_plan(&store, 1).unwrap();
        assert!(lean.frames.iter().all(|frame| frame.aa_seq.is_empty()));
        let subject = XPrelimSubject::Nt(SubjectNt::Packed(&store));
        for frame_index in 0..6 {
            let aa_len = full.frames[frame_index].aa_len;
            let buffer: Vec<u8> = full.frames[frame_index + 1..]
                .iter()
                .flat_map(|frame| frame.aa_seq[1..].iter().copied())
                .collect();
            for tail in 0..(HSP_MAX_WINDOW - 1).min(aa_len) {
                let chunk = TranslatedChunk::new(0, aa_len - tail);
                let after = |frames: &[QueryFrame]| {
                    XOneHitBuffers::new(frames, frame_index, chunk, subject, 1)
                        .unwrap()
                        .after_frame
                };
                let (from_lean, from_full) = (after(&lean.frames), after(&full.frames));
                let case = format!("{length} bases, frame index {frame_index}, tail {tail}");
                assert_eq!(from_lean, from_full, "{case}");
                assert!(
                    from_full.len() >= buffer.len().min(HSP_MAX_WINDOW - 1),
                    "{case}"
                );
                assert_eq!(from_full[..], buffer[..from_full.len()], "{case}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:603-605,635-657,703-712,
// 102-109 (cited at `x_add_one_hit_hsp`)
// Two one-hit HSPs of length -1 whose gapped starts lie 2 and 3 letters past a chunk of 30
// letters (`s_max` 31): B ends at (q10, s32) (W/W, then W/P), A at (q10, s33) with (q9, s32)
// (W/W twice, then C/P); their windows sum above 0 through the W/W letters right of the starts,
// which the left extension does not read. Both query ranges hold the query tree's middle (10)
// and share no end point, so both go to the subject tree, where they always go right: NCBI
// splits the node [31, 31] into [32, 31] and never ends. The search stops there.
#[test]
fn two_one_hit_hsps_past_their_chunk_stop_where_ncbi_splits_without_end() {
    const W: u8 = 20;
    const C: u8 = 3;
    const D: u8 = 4;
    const P: u8 = 14;
    let mut query = [C; 20];
    query[9..16].fill(W);
    let mut letters = vec![D; 60];
    letters[31] = P;
    letters[32..39].fill(W);
    let frame = one_hit_frame(1, &letters);
    let query_frames = [one_hit_frame(0, &query)];
    let frames: Vec<QueryFrame> = std::iter::once(frame.clone())
        .chain(
            [2, 3, -1, -2, -3]
                .into_iter()
                .map(|frame| one_hit_frame(frame, &[D; 60])),
        )
        .collect();
    let chunk = TranslatedChunk::new(0, 30);
    let buffers = XOneHitBuffers::new(
        &frames,
        0,
        chunk,
        super::super::x_subject_store::XPrelimSubject::Nt(SubjectNt::Letters(&[])),
        1,
    )
    .unwrap();
    let hit = |s_start: i32| InitHsp {
        frame: 1,
        chunk_offset: 0,
        q_seed: 5,
        s_seed: s_start as u32,
        q_start: 5,
        s_start,
        length: -1,
        score: 4,
    };
    let gapped = |hits: &[InitHsp]| {
        score_gapped_chunk::<_, true>(
            &[&query],
            &query_frames,
            &[0],
            20,
            &frame,
            chunk,
            &frame.aa_seq[1..31],
            &frame.aa_seq[1..],
            hits,
            11,
            1,
            38,
            &[1],
            ScoringMatrix::Blosum62,
            &mut GapAlignScratch::new(),
            super::super::score::MatrixScore,
            Some(&buffers),
        )
    };
    // Each alone is NCBI's HSP past the chunk.
    let one = |s_start: i32| gapped(&[hit(s_start)]).unwrap();
    let ends = |hsps: Vec<(usize, GappedHsp)>| {
        hsps.into_iter()
            .map(|(_, hsp)| (hsp.score, hsp.q_start, hsp.q_end, hsp.s_start, hsp.s_end))
            .collect::<Vec<_>>()
    };
    assert_eq!(ends(one(27)), vec![(11, 10, 10, 32, 32)]);
    assert_eq!(ends(one(28)), vec![(22, 9, 10, 32, 33)]);
    // Together they stop the search, in either order of the tree.
    let error = gapped(&[hit(27), hit(28)]).unwrap_err();
    assert!(
        error.to_string().contains("interval tree without end"),
        "{error}"
    );
    let mut a_first = hit(28);
    a_first.score = 5;
    let error = gapped(&[hit(27), a_first]).unwrap_err();
    assert!(
        error.to_string().contains("interval tree without end"),
        "{error}"
    );
}

// No NCBI counterpart: the pool of the units' buffers of a query batch (`XUnitResiduesPool`).
#[test]
fn x_unit_residues_pool_lends_its_idle_buffers_until_it_is_dropped() {
    let pool = XUnitResiduesPool::new();
    let mut first = pool.take();
    assert_eq!(
        first.1.capacity(),
        0,
        "an empty pool lends a new empty buffer"
    );
    first.1.resize(10, 7);
    let mut second = pool.take();
    assert_eq!(second.1.capacity(), 0, "a buffer in use is not lent");
    second.1.resize(20, 8);
    drop(first);
    drop(second);
    assert_eq!(pool.idle().len(), 2, "given back at the units' end");
    let unit = pool.take();
    assert_eq!(unit.1, vec![8; 20], "an idle buffer, as it was left");
    drop(unit);
    let (a, b) = (pool.take(), pool.take());
    let fresh = pool.take();
    assert_eq!(fresh.1.capacity(), 0, "both idle buffers are lent");
    drop(fresh);
    assert!(pool.idle().is_empty(), "an empty buffer is not kept");
    drop((a, b));
    assert_eq!(pool.idle().len(), 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2268-2323,2455-2537:
// s_QueryOffsetCompareHSPs returns zero when score and both endpoints
// match. The pinned NCBI Stage G probe retained the earlier score-355
// gapped_start (192, 47944) over the later (270, 48024).
#[test]
fn common_endpoint_purge_keeps_first_tied_gapped_seed() {
    let first = GappedHsp {
        frame: 3,
        score: 355,
        q_start: 107,
        q_end: 281,
        q_gapped_start: 192,
        s_start: 47852,
        s_end: 48034,
        s_gapped_start: 47944,
    };
    let later = GappedHsp {
        q_gapped_start: 270,
        s_gapped_start: 48024,
        ..first
    };
    let survivors = purge_preliminary_common_endpoints(vec![(0, first), (0, later)]);
    assert_eq!(survivors, vec![(0, first)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3924-3927:
// Read each query context's cutoff from the comparison-only probe.
fn multi_query_cutoffs() -> [i32; 3] {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/context_cutoffs_20260924/blosum62_word3.tsv"
    );
    let mut cutoffs = [None; 3];
    let mut count = 0;
    for line in fs::read_to_string(path).unwrap().lines() {
        let f: Vec<_> = line.split('\t').collect();
        let context: usize = f[2].parse().unwrap();
        let value: i32 = f[4].parse().unwrap();
        assert!(context < cutoffs.len());
        if let Some(previous) = cutoffs[context] {
            assert_eq!(previous, value);
        } else {
            cutoffs[context] = Some(value);
        }
        count += 1;
    }
    assert_eq!(count, 18);
    cutoffs.map(Option::unwrap)
}

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let text = fs::read_to_string(path).unwrap();
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(id) = line.strip_prefix('>') {
            records.push((id.to_string(), Vec::new()));
        } else {
            records.last_mut().unwrap().1.extend(line.as_bytes());
        }
    }
    records
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1810:
// compare the saved ordered NCBI link output after a Rust C-to-D call.
fn assert_stage_d_link_list_matches_trace(list: &LinkedHspList, trace: &str, event: usize) {
    let prefix = format!("D_HSP\t{event}\tlink_after\t");
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with(&prefix))
        .collect();
    assert_eq!(list.hsps.len(), expected.len());
    for (actual, line) in list.hsps.iter().zip(expected) {
        let f: Vec<_> = line.split('\t').collect();
        assert_eq!(actual.context, f[4].parse::<usize>().unwrap());
        assert_eq!(actual.hsp.frame, f[5].parse::<i8>().unwrap());
        assert_eq!(actual.hsp.q_start, f[6].parse::<i32>().unwrap());
        assert_eq!(actual.hsp.q_end, f[7].parse::<i32>().unwrap());
        assert_eq!(actual.hsp.s_start, f[8].parse::<i32>().unwrap());
        assert_eq!(actual.hsp.s_end, f[9].parse::<i32>().unwrap());
        assert_eq!(actual.hsp.score, f[10].parse::<i32>().unwrap());
        assert_eq!(actual.num, f[12].parse::<i32>().unwrap());
        assert_eq!(
            actual.evalue.to_bits(),
            f[13].parse::<f64>().unwrap().to_bits()
        );
    }
    let list_prefix = format!("D_LIST\t{event}\tlink_after\t");
    let list_row = trace
        .lines()
        .find(|line| line.starts_with(&list_prefix))
        .unwrap();
    assert_eq!(
        list.best_evalue.to_bits(),
        list_row
            .split('\t')
            .nth(6)
            .unwrap()
            .parse::<f64>()
            .unwrap()
            .to_bits()
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:964-985;
// c++/src/algo/blast/core/blast_parameters.c:774-815,902-999,302-383;
// c++/src/algo/blast/core/blast_engine.c:1390-1445,870-885:
// initial hit/word/link cutoff state is computed before the preliminary
// search; the complete appended list is passed to BLAST_LinkHsps.
#[test]
fn computed_default_parameters_preserve_ncbi_preliminary_link_input() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let gumbel = lookup_protein_gumbel_params(
        &ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        },
        (subject.len() / 3) as i64,
    )
    .unwrap();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:729-847:
    let parameters = local_parameters_for_call(
        &[(120, true), (70, true), (120, false)],
        subject.len() as i64,
        1,
        &[gapped; 3],
        &[ungapped; 3],
        LocalParameterOptions {
            expect_value: 10.0,
            do_sum_stats: true,
            max_intron_length: 0,
            gap_trigger_bits: 22.0,
            word_xdrop_bits: 7.0,
            scale_factor: 1.0,
            gumbel: Some(&gumbel),
            eff_searchsp_override: None,
            ungapped: None,
        },
        LocalParameterCall::Initial {
            min_subject_length: (subject.len() / 3) as i32,
            composition_based_stats: 2,
        },
    )
    .unwrap();
    let word_xdrop: Vec<_> = parameters.cutoffs.iter().map(|p| p.word_xdrop).collect();
    let word_cutoff: Vec<_> = parameters.cutoffs.iter().map(|p| p.word_cutoff).collect();
    let hit_cutoff: Vec<_> = parameters.cutoffs.iter().map(|p| p.hit_cutoff).collect();
    assert_eq!(hit_cutoff, [19, 17, i32::MAX]);
    assert_eq!(word_cutoff, [19, 17, i32::MAX]);
    assert_eq!(word_xdrop, [16, 16, 0]);
    let seg = SegParams::default();
    let profile = PreliminaryProfile {
        seg: Some(&seg),
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &word_xdrop,
        word_cutoff: &word_cutoff,
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &hit_cutoff,
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let preliminary = preliminary_hsps_of_letters(&refs, subject, 1, profile).unwrap();
    let trace = fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/run_20260924/multi_query_20260924_default.trace"
    ))
    .unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("D_HSP\t0\tlink_before\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[4].parse::<usize>().unwrap(),
                f[5].parse::<i8>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
                f[9].parse::<i32>().unwrap(),
                f[10].parse::<i32>().unwrap(),
            )
        })
        .collect();
    let actual: Vec<_> = preliminary
        .iter()
        .map(|(context, hsp)| {
            (
                *context,
                hsp.frame,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
                hsp.score,
            )
        })
        .collect();
    assert_eq!(actual, expected);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1810 and c++/src/algo/blast/core/blast_engine.c:643-676:
    // the scored Stage C list is linked before preliminary E-value reap.
    let mut linked = link_preliminary_hsps(
        &preliminary,
        &[120, 70, 120],
        &parameters.lengths,
        subject.len() as i32,
        &[gapped; 3],
        Some(&gumbel),
        &parameters.link.unwrap(),
    )
    .unwrap();
    assert_stage_d_link_list_matches_trace(&linked, &trace, 0);
    reap_by_evalue(&mut linked, parameters.prelim_evalue);
    assert_eq!(linked.hsps.len(), 20);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:964-985;
// c++/src/algo/blast/core/blast_engine.c:478-586,870-885:
// the initial computed hit/word cutoffs feed the chunked preliminary search;
// the appended output is the link input before NCBI modifies num/E-value.
#[test]
fn computed_positive_uneven_gap_inputs_match_ncbi_link_before() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/uneven_gap_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let gumbel = lookup_protein_gumbel_params(
        &ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        },
        (subject.len() / 3) as i64,
    )
    .unwrap();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:729-847:
    let parameters = local_parameters_for_call(
        &[(query.len(), true)],
        subject.len() as i64,
        1,
        &[gapped],
        &[ungapped],
        LocalParameterOptions {
            expect_value: 10.0,
            do_sum_stats: true,
            max_intron_length: 0,
            gap_trigger_bits: 22.0,
            word_xdrop_bits: 7.0,
            scale_factor: 1.0,
            gumbel: Some(&gumbel),
            eff_searchsp_override: None,
            ungapped: None,
        },
        LocalParameterCall::Initial {
            min_subject_length: (subject.len() / 3) as i32,
            composition_based_stats: 2,
        },
    )
    .unwrap();
    let word_xdrop = [parameters.cutoffs[0].word_xdrop];
    let word_cutoff = [parameters.cutoffs[0].word_cutoff];
    let gapped_cutoff = [parameters.cutoffs[0].hit_cutoff];
    assert_eq!(word_cutoff, [10]);
    let seg = SegParams::default();
    let profile = PreliminaryProfile {
        seg: Some(&seg),
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &word_xdrop,
        word_cutoff: &word_cutoff,
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &gapped_cutoff,
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let preliminary = preliminary_hsps_of_letters(&[query], subject, 1, profile).unwrap();
    let trace = fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/uneven_gap_run_20260924/uneven_gap_20260924_default.trace"
    ))
    .unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("D_HSP\t0\tlink_before\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[4].parse::<usize>().unwrap(),
                f[5].parse::<i8>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
                f[9].parse::<i32>().unwrap(),
                f[10].parse::<i32>().unwrap(),
            )
        })
        .collect();
    let actual: Vec<_> = preliminary
        .iter()
        .map(|(context, hsp)| {
            (
                *context,
                hsp.frame,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
                hsp.score,
            )
        })
        .collect();
    assert_eq!(actual, expected);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1810 and c++/src/algo/blast/core/blast_engine.c:643-676:
    // the scored Stage C list is linked before preliminary E-value reap.
    let mut linked = link_preliminary_hsps(
        &preliminary,
        &[query.len() as i32],
        &parameters.lengths,
        subject.len() as i32,
        &[gapped],
        Some(&gumbel),
        &parameters.link.unwrap(),
    )
    .unwrap();
    assert_stage_d_link_list_matches_trace(&linked, &trace, 0);
    reap_by_evalue(&mut linked, parameters.prelim_evalue);
    assert_eq!(linked.hsps.len(), 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2371-2468:
// s_AdjustHspOffsetsAndGetQueryData converts global initial offsets to
// local query coordinates and preserves the context on Blast_HSPInit.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3908-3919,4057-4091:
// the frame-local tree checks containment before each saved gapped HSP.
#[test]
fn multi_query_gapped_hsps_match_ncbi_context_scores_offsets_and_order() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let gapped_cutoff = multi_query_cutoffs();
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16, 16, 0],
        word_cutoff: &[0, 0, i32::MAX],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &gapped_cutoff,
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&refs, subject, 1, profile).unwrap();
    let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("GAPPED_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            assert_eq!(f[5], "0");
            (
                f[4].parse::<usize>().unwrap(),
                GappedHsp {
                    frame: f[9].parse().unwrap(),
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[10].parse().unwrap(),
                    s_end: f[11].parse().unwrap(),
                    s_gapped_start: f[12].parse().unwrap(),
                },
            )
        })
        .collect();
    let actual = stages
        .gapped
        .iter()
        .map(|(context, hsp, _)| (*context, *hsp))
        .collect::<Vec<_>>();
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(actual.len().min(expected.len()));
        panic!(
            "first multi-query gapped difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
            actual.get(first),
            expected.get(first),
            actual.len(),
            expected.len()
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:539-552,840-850:
// each frame is endpoint-purged, score-sorted, then appended to the
// cumulative list at kHspNumMax. Compare every append snapshot.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2809-2864:
// the first list is adopted and later lists are combined by score.
#[test]
fn multi_query_preliminary_append_snapshots_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let gapped_cutoff = multi_query_cutoffs();
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16, 16, 0],
        word_cutoff: &[0, 0, i32::MAX],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &gapped_cutoff,
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&refs, subject, 1, profile).unwrap();
    let actual = stages.appended.clone();
    let trace = fs::read_to_string(format!("{root}append_events.tsv")).unwrap();
    let mut expected: Vec<Vec<(usize, i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); 6];
    for line in trace
        .lines()
        .filter(|line| line.starts_with("APPEND_OUT_HSP\t"))
    {
        let f: Vec<_> = line.split('\t').collect();
        expected[f[1].parse::<usize>().unwrap()].push((
            f[3].parse().unwrap(),
            f[4].parse().unwrap(),
            f[5].parse().unwrap(),
            f[6].parse().unwrap(),
            f[7].parse().unwrap(),
            f[8].parse().unwrap(),
            f[9].parse().unwrap(),
        ));
    }
    let actual: Vec<Vec<_>> = actual
        .iter()
        .map(|snapshot| {
            snapshot
                .iter()
                .map(|(context, hsp)| {
                    (
                        *context,
                        hsp.score,
                        hsp.frame,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                    )
                })
                .collect()
        })
        .collect();
    assert_eq!(actual, expected);
    assert_eq!(
        actual.iter().map(Vec::len).collect::<Vec<_>>(),
        [3, 5, 11, 19, 25, 31]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:539-552,840-850:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2537,2809-2864:
// Remove a lower-scoring
// common-start HSP, then cap the cumulative score-sorted list.
#[test]
fn real_chunk_endpoint_deletion_and_append_cap_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let queries = read_fasta(&format!("{root}multi_query_20260924/query.faa"));
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let subject = &read_fasta(&format!("{root}multi_query_20260924/subjects.fna"))[0].1;
    let trace = fs::read_to_string(format!(
        "{root}chunk_purge_cap_real_path_20260924/trace.tsv"
    ))
    .unwrap();
    let mut expected: Vec<Vec<(usize, i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); 6];
    let mut first_input = Vec::new();
    let mut first_output = Vec::new();
    let mut purge_call = 0;
    for line in trace.lines() {
        let f: Vec<_> = line.split('\t').collect();
        match f[0] {
            "ENDPOINT_INPUT" => purge_call += 1,
            "ENDPOINT_IN_HSP" if purge_call == 1 => first_input.push((
                f[2].parse::<i32>().unwrap(),
                f[3].parse::<i8>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
            )),
            "ENDPOINT_OUT_HSP" if purge_call == 1 => first_output.push((
                f[2].parse::<i32>().unwrap(),
                f[3].parse::<i8>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
            )),
            "APPEND_OUT_HSP" => expected[f[1].parse::<usize>().unwrap()].push((
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
                f[9].parse().unwrap(),
            )),
            _ => {}
        }
    }
    assert_eq!(
        first_input,
        [(656, 1, 0, 200), (393, 1, 0, 225), (20, 1, 0, 200)]
    );
    assert_eq!(first_output, [(656, 1, 0, 200), (393, 1, 0, 225)]);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:525-552,840-850:
    // the comparison-only C probe changes this one HSP after the first
    // all six real append calls receive hsp_num_max=3.
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16, 16, 0],
        word_cutoff: &[0, 0, i32::MAX],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0, 0, i32::MAX],
        hsp_num_max: 3,
        gapped: true,
        subject_besthit: false,
    };
    let changes = std::rc::Rc::new(std::cell::Cell::new(0usize));
    let hook = {
        let changes = std::rc::Rc::clone(&changes);
        move |frame: i8, offset: usize, hsps: &mut Vec<(usize, GappedHsp)>| -> Result<()> {
            if frame == 1 && offset == 0 {
                let matches: Vec<_> = hsps
                    .iter_mut()
                    .filter(|(context, hsp)| {
                        *context == 0 && hsp.score == 20 && hsp.q_start == 17 && hsp.s_start == 1703
                    })
                    .collect();
                assert_eq!(matches.len(), 1);
                matches.into_iter().next().unwrap().1.q_start = 0;
                hsps.iter_mut()
                    .find(|(context, hsp)| {
                        *context == 0 && hsp.score == 20 && hsp.q_start == 0 && hsp.s_start == 1703
                    })
                    .unwrap()
                    .1
                    .s_start = 200;
                changes.set(changes.get() + 1);
            }
            Ok(())
        }
    };
    let (_, events) =
        preliminary_hsps_traced(&refs, subject, 1, profile, Some(Box::new(hook))).unwrap();
    assert_eq!(changes.get(), 1);
    let first_purged = events
        .iter()
        .find_map(|event| match event {
            PreliminaryEvent::Purged(1, 0, hsps) => Some(hsps),
            _ => None,
        })
        .unwrap();
    let actual_purged: Vec<_> = first_purged
        .iter()
        .map(|(_, hsp)| (hsp.score, hsp.frame, hsp.q_start, hsp.s_start))
        .collect();
    assert_eq!(actual_purged, first_output);
    let mut appended = Vec::new();
    for event in &events {
        if let PreliminaryEvent::Appended(_, hsps) = event {
            appended.push(
                hsps.iter()
                    .map(|(context, hsp)| {
                        (
                            *context,
                            hsp.score,
                            hsp.frame,
                            hsp.q_start,
                            hsp.q_end,
                            hsp.s_start,
                            hsp.s_end,
                        )
                    })
                    .collect::<Vec<_>>(),
            );
        }
    }
    assert_eq!(appended, expected);
    assert_eq!(
        appended.iter().map(Vec::len).collect::<Vec<_>>(),
        [2, 3, 3, 3, 3, 3]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:259-312,635-721:
// query-indexed preliminary lists enter traceback independently;
// a fence retries the same original list on the full subject.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1644-1684:
// only the second return with fence_hit=0 is retained here.
#[test]
fn multi_query_traceback_lists_match_ncbi_post_fence_hsps() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let gapped_cutoff = multi_query_cutoffs();
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16, 16, 0],
        word_cutoff: &[0, 0, i32::MAX],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &gapped_cutoff,
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&refs, subject, 1, profile).unwrap();
    let appended = stages.list.clone();

    let trace = fs::read_to_string(format!("{root}traceback_events.tsv")).unwrap();
    let mut expected: Vec<Vec<(i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); refs.len()];
    let mut context = 0usize;
    let mut successful_pass = false;
    for line in trace.lines() {
        let f: Vec<_> = line.split('\t').collect();
        match f[0] {
            "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
            "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
            "TRACEBACK_OUT_HSP" if successful_pass => expected[context].push((
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
            )),
            _ => {}
        }
    }
    for context in [1usize, 0] {
        let per_query: Vec<_> = appended
            .iter()
            .filter(|(index, _)| *index == context)
            .map(|(_, hsp)| *hsp)
            .collect();
        let actual =
            full_translation_traceback_blosum62(refs[context], subject, 1, &per_query, 11, 1, 64)
                .unwrap();
        let actual: Vec<_> = actual
            .into_iter()
            .map(|(hsp, full_retry)| {
                assert!(full_retry, "NCBI retried this list after a fence");
                (
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
            .collect();
        assert_eq!(actual, expected[context], "query context {context}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-605:
// The comparison-only C API input sets percent_identity=100.0 during
// each real NCBI HSPTest call, yielding positive deletions.
#[test]
fn real_path_hsp_test_deletion_order_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let queries = read_fasta(&format!("{root}multi_query_20260924/query.faa"));
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let subject = &read_fasta(&format!("{root}multi_query_20260924/subjects.fna"))[0].1;
    let gapped_cutoff = multi_query_cutoffs();
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16, 16, 0],
        word_cutoff: &[0, 0, i32::MAX],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &gapped_cutoff,
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&refs, subject, 1, profile).unwrap();
    let appended = stages.list.clone();
    let trace = fs::read_to_string(format!("{root}hsp_test_real_path_20260924/trace.tsv")).unwrap();
    let mut expected_tests = vec![Vec::new(); refs.len()];
    let mut expected_output = vec![Vec::new(); refs.len()];
    let mut context = 0usize;
    let mut successful_pass = false;
    for line in trace.lines() {
        let f: Vec<_> = line.split('\t').collect();
        match f[0] {
            "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
            "HSP_TEST" => expected_tests[context].push((
                f[1] == "1",
                f[2].parse::<usize>().unwrap(),
                f[3].parse::<usize>().unwrap(),
                f[4].parse::<usize>().unwrap(),
                f[5].parse::<usize>().unwrap(),
                f[6].parse::<usize>().unwrap(),
            )),
            "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
            "TRACEBACK_OUT_HSP" if successful_pass => expected_output[context].push((
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i8>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            )),
            _ => {}
        }
    }
    assert_eq!(expected_tests.iter().map(Vec::len).sum::<usize>(), 31);
    assert_eq!(
        expected_tests.iter().flatten().filter(|row| row.0).count(),
        17
    );
    for context in [1usize, 0] {
        let per_query: Vec<_> = appended
            .iter()
            .filter(|(index, _)| *index == context)
            .map(|(_, hsp)| *hsp)
            .collect();
        let mut events = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events(
            refs[context],
            subject,
            1,
            &per_query,
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            100.0,
            None,
            Some(&mut events),
            None,
            None,
            None,
        )
        .unwrap();
        assert_eq!(
            events, expected_tests[context],
            "query context {context} HSPTest order"
        );
        let actual: Vec<_> = actual
            .into_iter()
            .map(|(hsp, retried)| {
                assert!(retried);
                (
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
            .collect();
        assert_eq!(
            actual, expected_output[context],
            "query context {context} survivors"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3259-3321:
// q_length <= 11 returns the midpoint; a strictly positive sliding
// maximum wins; otherwise the terminal 11-residue window is tested.
#[test]
fn protein_start_offset_window_boundaries_follow_ncbi() {
    let encode = |residues: &[u8]| {
        let frame = encode_protein_query_frame_with_seg(residues, None);
        frame.aa_seq[1..1 + residues.len()].to_vec()
    };
    let all_m = encode(&b"M".repeat(20));
    assert_eq!(
        blast_get_offsets_for_gapped_alignment_protein(
            &all_m[..11],
            &all_m[..11],
            0,
            11,
            0,
            11,
            ScoringMatrix::Blosum62,
        ),
        Some((5, 5))
    );
    assert_eq!(
        blast_get_offsets_for_gapped_alignment_protein(
            &all_m,
            &all_m,
            0,
            20,
            0,
            20,
            ScoringMatrix::Blosum62,
        ),
        Some((10, 10))
    );
    let mut end_only = b"D".repeat(9);
    end_only.extend_from_slice(&b"M".repeat(11));
    let end_only = encode(&end_only);
    assert_eq!(
        blast_get_offsets_for_gapped_alignment_protein(
            &end_only,
            &all_m[..12],
            0,
            20,
            0,
            12,
            ScoringMatrix::Blosum62,
        ),
        Some((15, 7))
    );
    let all_k = encode(&b"K".repeat(20));
    assert_eq!(
        blast_get_offsets_for_gapped_alignment_protein(
            &all_k,
            &all_m,
            0,
            20,
            0,
            20,
            ScoringMatrix::Blosum62,
        ),
        None
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-445:
// Recomputes starting offsets via `BlastGetOffsetsForGappedAlignment`; when that fails,
// frees the HSP (`Blast_HSPFree`) and continues to the next entry.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3259-3321: positive 11-residue window scores
// choose a start; two nonpositive endpoint windows return FALSE.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:446-448: on TRUE, write q_start and s_start
// into the HSP before alignment; the real-path positive probe tests this.
fn compare_real_path_start_offset_with_ncbi(positive: bool) {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let queries = read_fasta(&format!("{root}multi_query_20260924/query.faa"));
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let subject = &read_fasta(&format!("{root}multi_query_20260924/subjects.fna"))[0].1;
    let gapped_cutoff = multi_query_cutoffs();
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16, 16, 0],
        word_cutoff: &[0, 0, i32::MAX],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &gapped_cutoff,
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&refs, subject, 1, profile).unwrap();
    let appended = stages.list.clone();
    let case = if positive {
        "start_success_real_path_20260924"
    } else {
        "start_failure_real_path_20260924"
    };
    let trace = fs::read_to_string(format!("{root}{case}/trace.tsv")).unwrap();
    let mut expected_output = vec![Vec::new(); refs.len()];
    let mut expected_starts = Vec::new();
    let mut expected_after = Vec::new();
    let mut context = 0usize;
    let mut successful_pass = false;
    for line in trace.lines() {
        let f: Vec<_> = line.split('\t').collect();
        match f[0] {
            "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
            "START_RESULT" => expected_starts.push((
                f[1] == "1",
                f[2].parse::<i32>().unwrap(),
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
            )),
            "AFTER_START_HSP" => expected_after.push((
                f[2].parse::<i32>().unwrap(),
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            )),
            "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
            "TRACEBACK_OUT_HSP" if successful_pass => expected_output[context].push((
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i8>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            )),
            _ => {}
        }
    }
    assert_eq!(expected_starts.len(), 2);
    assert!(expected_starts.iter().all(|row| row.0 == positive));
    for context in [1usize, 0] {
        let mut per_query: Vec<_> = appended
            .iter()
            .filter(|(index, _)| *index == context)
            .map(|(_, hsp)| *hsp)
            .collect();
        if context == 0 {
            // Comparison-only NCBI input-state intervention in
            // ncbi_start_failure_inject.c, before the traceback call.
            per_query[0].q_start = 0;
            per_query[0].q_end = if positive { 120 } else { 20 };
            per_query[0].q_gapped_start = 0;
            per_query[0].frame = 1;
            per_query[0].s_start = if positive { 200 } else { 0 };
            per_query[0].s_end = if positive { 320 } else { 20 };
            per_query[0].s_gapped_start = 0;
        }
        let mut starts = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events(
            refs[context],
            subject,
            1,
            &per_query,
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            0.0,
            None,
            None,
            Some(&mut starts),
            None,
            None,
        )
        .unwrap();
        if context == 0 {
            assert_eq!(starts, expected_starts);
        } else {
            assert!(starts.is_empty());
        }
        if context == 0 && positive {
            let hsp = &actual[0].0;
            let &(score, q_start, q_end, q_gap, s_start, s_end, s_gap) =
                expected_after.last().unwrap();
            assert_eq!(
                (
                    hsp.score,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.q_gapped_start,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.s_gapped_start
                ),
                (score, q_start, q_end, q_gap, s_start, s_end, s_gap),
                "NCBI successful acquired start writeback"
            );
        }
        let actual: Vec<_> = actual
            .into_iter()
            .map(|(hsp, retried)| {
                assert!(retried);
                (
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
            .collect();
        assert_eq!(actual, expected_output[context], "query context {context}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-445:
// When start-offset lookup fails, frees the HSP (`Blast_HSPFree`) and skips it.
#[test]
fn real_path_start_offset_failure_deletion_matches_ncbi() {
    compare_real_path_start_offset_with_ncbi(false);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-448:
#[test]
fn real_path_start_offset_success_writeback_matches_ncbi() {
    compare_real_path_start_offset_with_ncbi(true);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:572-586:
// adjust the new HSP's subject offsets, then merge with overlap 100.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-3035:
// only equal-context, close-diagonal overlap HSPs are merged.
// The pinned oracle records old/new lists immediately before each call.
#[test]
fn long_chunk_overlap_merge_matches_ncbi_input_score_and_order() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/long_chunk_20260924/"
    );
    let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let parse_hsp = |line: &str| {
        let f: Vec<_> = line.split('\t').collect();
        (
            f[3].parse::<usize>().unwrap(),
            GappedHsp {
                score: f[4].parse().unwrap(),
                frame: f[5].parse().unwrap(),
                q_start: f[6].parse().unwrap(),
                q_end: f[7].parse().unwrap(),
                q_gapped_start: f[8].parse().unwrap(),
                s_start: f[9].parse().unwrap(),
                s_end: f[10].parse().unwrap(),
                s_gapped_start: f[11].parse().unwrap(),
            },
        )
    };
    let rows: Vec<_> = trace.lines().collect();
    let for_call = |event: &str, call: usize| {
        rows.iter()
            .copied()
            .filter(|line| {
                let mut fields = line.split('\t');
                fields.next() == Some(event)
                    && fields.next() == Some(if call == 0 { "0" } else { "1" })
            })
            .map(&parse_hsp)
            .collect::<Vec<_>>()
    };
    let first = merge_subject_chunk_hsps(
        Vec::new(),
        for_call("MERGE_IN_HSP", 0),
        0,
        0,
        i32::MAX as usize,
        true,
    );
    assert_eq!(first, for_call("MERGE_OUT_HSP", 0));
    let old = for_call("MERGE_OLD_HSP", 1);
    assert_eq!(first, old);
    let second = merge_subject_chunk_hsps(
        old,
        for_call("MERGE_IN_HSP", 1),
        4_999_900,
        100,
        i32::MAX as usize,
        true,
    );
    assert_eq!(second, for_call("MERGE_OUT_HSP", 1));
    assert_eq!(second[0].1.score, 660);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2762-2864:
// when its five compared fields tie, choose the old list's HSP first.
// The expected rows come from the direct comparison-only NCBI C API.
#[test]
fn preliminary_append_cap_and_tie_match_ncbi_api() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/append_cap_20260924/oracle.tsv"
    );
    let hsp = |context: usize, score: i32, frame: i8, start: i32| {
        (
            context,
            GappedHsp {
                frame,
                score,
                q_start: 0,
                q_end: 20,
                q_gapped_start: 5,
                s_start: start,
                s_end: start + 20,
                s_gapped_start: start + 5,
            },
        )
    };
    let input = [
        hsp(0, 100, 1, 100),
        hsp(0, 80, 1, 300),
        hsp(0, 60, 1, 500),
        hsp(1, 90, 2, 200),
        hsp(1, 80, 2, 300),
        hsp(1, 70, 2, 400),
    ];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:539-552,840-850: each frame's list is
    // endpoint-purged, then appended to the cumulative list at hsp_num_max = 3.
    let of_frame = |frame: i8| -> Vec<_> {
        input
            .iter()
            .copied()
            .filter(|(_, hsp)| hsp.frame == frame)
            .collect()
    };
    let first = append_preliminary_hsps(
        Vec::new(),
        purge_preliminary_common_endpoints(of_frame(1)),
        3,
    );
    let snapshot =
        append_preliminary_hsps(first, purge_preliminary_common_endpoints(of_frame(2)), 3);
    let actual: Vec<_> = snapshot
        .iter()
        .enumerate()
        .map(|(index, (context, hsp))| {
            (
                index,
                *context,
                hsp.score,
                hsp.frame,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
            )
        })
        .collect();
    let trace = fs::read_to_string(path).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("CAP_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].parse::<usize>().unwrap(),
                f[2].parse::<usize>().unwrap(),
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i8>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            )
        })
        .collect();
    assert_eq!(actual, expected);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3936-3953,4071-4091:
// Finds the gapped start via `BlastGetStartForGappedAlignment` and runs
// `s_BlastProtGappedAlignment`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:811-811
fn compare_fixture(fixture_name: &str, oracle_name: &str, lowercase: bool) {
    let root = format!(
        "{}/../docs/evidence/tlosan_stage_c/{fixture_name}/",
        env!("CARGO_MANIFEST_DIR")
    );
    let oracle = format!(
        "{}/../docs/evidence/tlosan_stage_c/{oracle_name}/gapped_output.tsv",
        env!("CARGO_MANIFEST_DIR")
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subjects = read_fasta(&format!("{root}subjects.fna"));
    let mut actual = Vec::new();
    for (id, subject) in &subjects {
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: lowercase,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
            gapped: true,
            subject_besthit: false,
        };
        let stages = preliminary_stages(&[query], subject, 1, profile).unwrap();
        for (_, hsp, _) in stages.gapped {
            actual.push((id.clone(), hsp));
        }
    }
    let expected: Vec<_> = fs::read_to_string(oracle)
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].to_string(),
                GappedHsp {
                    frame: f[2].parse().unwrap(),
                    score: f[4].parse().unwrap(),
                    q_start: f[7].parse().unwrap(),
                    q_end: f[8].parse().unwrap(),
                    q_gapped_start: f[9].parse().unwrap(),
                    s_start: f[11].parse().unwrap(),
                    s_end: f[12].parse().unwrap(),
                    s_gapped_start: f[13].parse().unwrap(),
                },
            )
        })
        .collect();
    assert_eq!(actual, expected);
}

#[test]
fn fixed_local_subject_gapped_scores_offsets_and_order_match_ncbi() {
    compare_fixture("run_20260923", "gapped_20260924", false);
}

#[test]
fn lowercase_gapped_scores_offsets_and_order_match_ncbi() {
    compare_fixture("lowercase_20260923", "gapped_lowercase_20260924", true);
}

#[test]
fn ambiguity_gapped_scores_offsets_and_order_match_ncbi() {
    compare_fixture("ambiguity_20260923", "gapped_ambiguity_20260924", false);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-552:
// Runs `GetGappedScore` on the initial hit list to build the gapped HSP list.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3936-3953:
// Finds the gapped start via `BlastGetStartForGappedAlignment` and runs
// `s_BlastProtGappedAlignment`.
#[test]
fn code32_api_gapped_scores_offsets_and_order_match_ncbi() {
    let fixture = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_a/fixtures/"
    );
    let oracle = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/gapped_code32_20260924/gapped_code32_api.tsv"
    );
    let query = &read_fasta(&format!("{fixture}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{fixture}subject_code32.fna"))[0].1;
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[13],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[1],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&[query], subject, 32, profile).unwrap();
    let actual = stages
        .gapped
        .iter()
        .map(|(_, hsp, _)| *hsp)
        .collect::<Vec<_>>();
    let expected: Vec<GappedHsp> = fs::read_to_string(oracle)
        .unwrap()
        .lines()
        .skip(1)
        .filter_map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            if f[0] != "32" {
                return None;
            }
            Some(GappedHsp {
                frame: f[16].parse().unwrap(),
                score: f[10].parse().unwrap(),
                q_start: f[13].parse().unwrap(),
                q_end: f[14].parse().unwrap(),
                q_gapped_start: f[15].parse().unwrap(),
                s_start: f[17].parse().unwrap(),
                s_end: f[18].parse().unwrap(),
                s_gapped_start: f[19].parse().unwrap(),
            })
        })
        .collect();
    assert_eq!(actual, expected);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
    // same function order and code-32 subject translation enter the
    // comparison-only NCBI C++ API gapped score oracle.
    assert_eq!(
        stages.list,
        expected.into_iter().map(|hsp| (0, hsp)).collect::<Vec<_>>()
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1644-1677:
// Refetches the whole subject when traceback hits a fence sentinel.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-589:
#[test]
fn full_translation_traceback_matches_ncbi_post_fence_scores_and_offsets() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/run_20260923/"
    );
    let oracle = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/gapped_20260924/traceback_hsps.tsv"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subjects = read_fasta(&format!("{root}subjects.fna"));
    let mut actual = Vec::new();
    for (id, subject) in &subjects {
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
            gapped: true,
            subject_besthit: false,
        };
        let stages = preliminary_stages(&[query], subject, 1, profile).unwrap();
        let gapped = stages
            .gapped
            .iter()
            .map(|(_, hsp, _)| *hsp)
            .collect::<Vec<_>>();
        for (hsp, fence_hit) in
            full_translation_traceback_blosum62(query, subject, 1, &gapped, 11, 1, 64).unwrap()
        {
            actual.push((
                id.clone(),
                hsp.frame,
                hsp.score,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
                fence_hit,
            ));
        }
    }
    let expected: Vec<_> = fs::read_to_string(oracle)
        .unwrap()
        .lines()
        .skip(1)
        .filter_map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            if f[3] != "2" {
                return None;
            }
            Some((
                f[2].to_string(),
                f[16].parse().unwrap(),
                f[15].parse().unwrap(),
                f[17].parse().unwrap(),
                f[18].parse().unwrap(),
                f[19].parse().unwrap(),
                f[20].parse().unwrap(),
                true,
            ))
        })
        .collect();
    assert_eq!(actual, expected);
    assert!(actual
        .iter()
        .any(|(id, _, score, _, _, _, _, _)| { id == "ambiguous" && *score == 646 }));
}

// No NCBI counterpart: the translated windows of a subject held packed (EXPERIMENT,
// LOSAT_X_TBNSTORE) are those of its letters: partial windows of both strands, one read again,
// one read by the other forward frames, windows that overlap, a window at the end of the subject and the whole subject (after a
// fence) of a forward and a reverse frame, under genetic codes 1 and 4.
#[test]
fn target_translation_of_a_packed_subject_is_that_of_its_letters() {
    let mut state = 11u64;
    let letters: Vec<u8> = (0..3_000)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let draw = (state >> 33) as usize;
            if draw.is_multiple_of(16) {
                b"RYKMSWBDHVN-rykmswbdhvn"[draw / 16 % 23]
            } else {
                b"ACGTacgt"[draw / 16 % 8]
            }
        })
        .collect();
    let store = crate::sequence::packed_subject::PackedNucleotide::from_subject_letters(&letters)
        .expect("letters with ncbi4na codes");
    let requests: [(i8, i32, i32); 12] = [
        (1, 10, 60),
        (2, 10, 60),
        (3, 10, 60),
        (-2, 500, 560),
        (-2, 500, 560),
        (1, 700, 760),
        (1, 720, 790),
        (3, 0, 30),
        (-1, 950, 990),
        (2, -1, 0),
        (-3, -1, 0),
        (2, 20, 40),
    ];
    for id in [1, 4] {
        let code = GeneticCode::try_from_id(id).unwrap();
        let mut from_letters = TargetTranslation::of(SubjectNt::Letters(&letters), &code);
        let mut from_store = TargetTranslation::of(SubjectNt::Packed(&store), &code);
        for (frame, offset, end) in requests {
            let expected = from_letters
                .get(frame, offset, end)
                .map(|(sequence, stop, base)| (sequence.to_vec(), stop, base))
                .unwrap();
            let got = from_store
                .get(frame, offset, end)
                .map(|(sequence, stop, base)| (sequence.to_vec(), stop, base))
                .unwrap();
            assert_eq!(got, expected, "code {id}, frame {frame}, {offset}..{end}");
        }
    }
}

// No NCBI counterpart: the traceback of a subject held packed (EXPERIMENT, LOSAT_X_TBNSTORE)
// reads windows of its store, made in upper case, where the record's letters were read. The
// subjects of the post-fence run (each list retried with the whole subject of its frames, both
// strands) and of the long multi-HSP run (partial windows), in lower case with ambiguous
// letters and `-`, under genetic codes 1 and 4.
#[test]
fn full_translation_traceback_reads_a_packed_subject_as_its_letters() {
    let evidence = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let (mut retried, mut minus, mut rows) = (false, false, 0);
    for run in ["run_20260923", "multi_hsp_20260924"] {
        let query = &read_fasta(&format!("{evidence}{run}/query.faa"))[0].1;
        for (id, subject) in read_fasta(&format!("{evidence}{run}/subjects.fna")) {
            let mut letters = subject.to_ascii_lowercase();
            for (index, letter) in letters.iter_mut().enumerate() {
                if index % 37 == 11 {
                    *letter = b"RYKMSWBDHVN-"[index % 12];
                }
            }
            let store =
                crate::sequence::packed_subject::PackedNucleotide::from_subject_letters(&letters)
                    .expect("letters with ncbi4na codes");
            for db_gencode in [1, 4] {
                let profile = PreliminaryProfile {
                    seg: None,
                    soft_masking: false,
                    threshold: 13,
                    window: 40,
                    word_xdrop: &[16],
                    word_cutoff: &[0],
                    mask_lowercase: false,
                    matrix: ScoringMatrix::Blosum62,
                    word_size: 3,
                    gap_open: 11,
                    gap_extend: 1,
                    gap_xdrop: 38,
                    gapped_cutoff: &[0],
                    hsp_num_max: i32::MAX as usize,
                    gapped: true,
                    subject_besthit: false,
                };
                let stages = preliminary_stages(&[query], &letters, db_gencode, profile).unwrap();
                let mut gapped = stages
                    .gapped
                    .iter()
                    .map(|(_, hsp, _)| *hsp)
                    .collect::<Vec<_>>();
                sort_gapped_score_if_needed(&mut gapped);
                let traceback = |subject| {
                    full_translation_traceback_with_matrix_and_events_with_mask_mode_owned(
                        query,
                        subject,
                        db_gencode,
                        &gapped,
                        ScoringMatrix::Blosum62,
                        11,
                        1,
                        64,
                        0.0,
                        None,
                        false,
                        false,
                        None,
                        None,
                        None,
                        None,
                    )
                    .unwrap()
                };
                let from_letters = traceback(SubjectNt::Letters(&letters));
                let from_store = traceback(SubjectNt::Packed(&store));
                assert_eq!(
                    format!("{from_store:?}"),
                    format!("{from_letters:?}"),
                    "{run} {id}, genetic code {db_gencode}"
                );
                rows += from_letters.0.len();
                retried |= from_letters.0.iter().any(|row| row.retried);
                minus |= from_letters.0.iter().any(|row| row.hsp.frame < 0);
            }
        }
    }
    assert!(
        rows > 0 && retried && minus,
        "{rows} rows, retried {retried}, minus {minus}"
    );
}

// No NCBI counterpart: the windows of the whole translation that the traceback reads after a
// fence (`x_whole_window`) hold the bytes of the whole translation (`get` with offset -1)
// between their ends, end as it does at its ends and with a fence elsewhere, and give its
// length; windows of every frame, at the subject's ends, read again, inside a wider one and
// wider than the frame, of the letters and of the store, under genetic codes 1 and 4; a subject
// of fewer than 12 letters gives what `get(frame, -1, _)` gives it (no residues).
#[test]
fn whole_windows_hold_the_whole_translation() {
    let mut state = 23u64;
    let letters: Vec<u8> = (0..9_001)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let draw = (state >> 33) as usize;
            if draw.is_multiple_of(16) {
                b"RYKMSWBDHVN-rykmswbdhvn"[draw / 16 % 23]
            } else {
                b"ACGTacgt"[draw / 16 % 8]
            }
        })
        .collect();
    let store = crate::sequence::packed_subject::PackedNucleotide::from_subject_letters(&letters)
        .expect("letters with ncbi4na codes");
    let requests: [(i32, i32, usize); 9] = [
        (1_000, 1_100, HSP_MAX_WINDOW),
        (1_000, 1_100, HSP_MAX_WINDOW),
        (1_005, 1_090, 12),
        (0, 40, 13),
        (2_950, 2_985, HSP_MAX_WINDOW),
        (2_980, 2_999, 64),
        (1_500, 1_600, 300),
        (30, 2_950, HSP_MAX_WINDOW),
        (1_200, 1_300, usize::MAX),
    ];
    for id in [1, 4] {
        let code = GeneticCode::try_from_id(id).unwrap();
        for subject in [SubjectNt::Letters(&letters), SubjectNt::Packed(&store)] {
            let mut whole = TargetTranslation::of(subject, &code);
            let mut windows = TargetTranslation::of(subject, &code);
            for frame in [1, 2, 3, -1, -2, -3] {
                let (residues, length, base) = whole
                    .get(frame, -1, 0)
                    .map(|(sequence, stop, base)| (sequence.to_vec(), stop, base))
                    .unwrap();
                assert_eq!((residues.len(), base), (length + 1, 0));
                for (offset, end, margin) in requests {
                    let (sequence, got_length, got_base) = windows
                        .x_whole_window(frame, offset, end, margin)
                        .map(|(sequence, stop, base)| (sequence.to_vec(), stop, base))
                        .unwrap();
                    let what = format!("code {id} frame {frame} {offset}..{end} +-{margin}");
                    assert_eq!(got_length, length, "{what}");
                    // A fence (201, never a residue) before the window's first residue
                    // (subject position `got_base`), or the first residue of the frame.
                    let left = usize::from(sequence[0] == 201);
                    let first = if left == 1 {
                        usize::try_from(got_base + 1).unwrap()
                    } else {
                        assert_eq!(got_base, 0, "{what}");
                        0
                    };
                    let last = first + sequence.len() - left - 1;
                    assert!(
                        first <= (offset as usize).saturating_sub(margin)
                            && last >= (end as usize).saturating_add(margin).min(length),
                        "{what}: [{first}, {last})"
                    );
                    assert_eq!(
                        sequence[left..sequence.len() - 1],
                        residues[first..last],
                        "{what}"
                    );
                    let right = if last == length { 0 } else { 201 };
                    assert_eq!(*sequence.last().unwrap(), right, "{what}");
                }
            }
        }
    }
    let short = b"ATGGCCAAATT".to_vec();
    let code = GeneticCode::try_from_id(1).unwrap();
    for frame in [1, 2, 3, -1, -2, -3] {
        let mut whole = TargetTranslation::new(&short, &code);
        let mut windows = TargetTranslation::new(&short, &code);
        let expected = whole.get(frame, -1, 2).map(|_| ()).unwrap_err().to_string();
        let got = windows
            .x_whole_window(frame, 0, 2, HSP_MAX_WINDOW)
            .map(|_| ())
            .unwrap_err()
            .to_string();
        assert_eq!(got, expected, "frame {frame}");
        assert_eq!(
            got,
            x_traceback_of_an_untranslated_frame().to_string(),
            "frame {frame}"
        );
        // A partial window of the HSP refuses too (NCBI's `subject.offset >= 0` takes the same
        // skip).
        let partial = windows
            .get(frame, 0, 2)
            .map(|_| ())
            .unwrap_err()
            .to_string();
        assert_eq!(partial, expected, "frame {frame}");
    }
    // The bound: NCBI translates a subject of 12 bases (`0 < 12 / 3 - 3`), none shorter.
    for n in 1..=12 {
        let subject = vec![b'A'; n];
        for frame in [1, 2, 3, -1, -2, -3] {
            let read = TargetTranslation::new(&subject, &code)
                .get(frame, 0, 1)
                .is_ok();
            assert_eq!(read, n == 12, "{n} bases, frame {frame}");
        }
    }
}

// No NCBI counterpart: after a fence the traceback reads windows of the whole translation
// (`x_whole_window`), traced again with wider windows when one is too narrow. Its HSPs,
// `stat_length` and events are those of the whole translation (`usize::MAX`: every window is
// the whole frame, the layout of `get(frame, -1, _)`, `whole_windows_hold_the_whole_translation`)
// whatever the first windows: from `HSP_MAX_WINDOW` residues on either side of an HSP (a list
// retried after a fence is traced again in wider windows) to the default. The subjects of the
// post-fence and multi-HSP runs, and the same between long random flanks (where the default
// windows end inside the frame), letters and store, under genetic codes 1 and 4.
#[test]
fn whole_windows_trace_as_the_whole_translation() {
    let evidence = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let mut state = 31u64;
    let mut flank = |length: usize| -> Vec<u8> {
        (0..length)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 33) as usize % 4]
            })
            .collect()
    };
    let (mut retried, mut rows, mut flanked, mut widened) = (false, 0, false, false);
    for run in ["run_20260923", "multi_hsp_20260924"] {
        let query = &read_fasta(&format!("{evidence}{run}/query.faa"))[0].1;
        for (id, subject) in read_fasta(&format!("{evidence}{run}/subjects.fna")) {
            let mut padded = flank(9_001);
            padded.extend_from_slice(&subject);
            padded.extend(flank(8_999));
            for (letters, long) in [(subject.clone(), false), (padded, true)] {
                let store =
                    crate::sequence::packed_subject::PackedNucleotide::from_subject_letters(
                        &letters,
                    )
                    .expect("letters with ncbi4na codes");
                for db_gencode in [1, 4] {
                    let profile = PreliminaryProfile {
                        seg: None,
                        soft_masking: false,
                        threshold: 13,
                        window: 40,
                        word_xdrop: &[16],
                        word_cutoff: &[0],
                        mask_lowercase: false,
                        matrix: ScoringMatrix::Blosum62,
                        word_size: 3,
                        gap_open: 11,
                        gap_extend: 1,
                        gap_xdrop: 38,
                        gapped_cutoff: &[0],
                        hsp_num_max: i32::MAX as usize,
                        gapped: true,
                        subject_besthit: false,
                    };
                    let stages =
                        preliminary_stages(&[query], &letters, db_gencode, profile).unwrap();
                    let mut gapped = stages
                        .gapped
                        .iter()
                        .map(|(_, hsp, _)| *hsp)
                        .collect::<Vec<_>>();
                    sort_gapped_score_if_needed(&mut gapped);
                    let traceback = |subject, margin| {
                        let mut events = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
                        let (rows, stat_length, attempts) = x_traceback_in_windows(
                            query,
                            subject,
                            db_gencode,
                            &gapped,
                            ScoringMatrix::Blosum62,
                            11,
                            1,
                            64,
                            0.0,
                            None,
                            false,
                            false,
                            Some(&mut events.0),
                            Some(&mut events.1),
                            Some(&mut events.2),
                            Some(&mut events.3),
                            margin,
                            crate::algorithm::tblastn::score::MatrixScore,
                        )
                        .unwrap();
                        (format!("{rows:?} {stat_length} {events:?}"), rows, attempts)
                    };
                    let (expected, whole, once) =
                        traceback(SubjectNt::Letters(&letters), usize::MAX);
                    // whole frames are traced once after a fence
                    assert!(once <= 1, "{run} {id}: {once} attempts in whole frames");
                    for subject in [SubjectNt::Letters(&letters), SubjectNt::Packed(&store)] {
                        for margin in [HSP_MAX_WINDOW, 12, 33, 300, X_WHOLE_MARGIN, usize::MAX] {
                            let (got, _, attempts) = traceback(subject, margin);
                            assert_eq!(
                                got, expected,
                                "{run} {id} long {long}, code {db_gencode}, margin {margin}"
                            );
                            assert_eq!(attempts > 0, once > 0, "{run} {id} margin {margin}");
                            widened |= attempts > 1;
                        }
                    }
                    rows += whole.len();
                    if whole.iter().any(|row| row.retried) {
                        retried = true;
                        flanked |= long;
                    }
                }
            }
        }
    }
    assert!(
        rows > 0 && retried && flanked && widened,
        "{rows} rows, retried {retried}, retried with flanks {flanked}, windows widened {widened}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1154-1228:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1644-1684:
// Refetches the whole subject and retries the HSP list when traceback hits a fence sentinel.
#[test]
fn long_multi_hsp_partial_windows_and_retry_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_hsp_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let code = GeneticCode::try_from_id(1).unwrap();
    let mut target = TargetTranslation::new(subject, &code);
    for (index, line) in fs::read_to_string(format!("{root}translation_windows.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .enumerate()
    {
        if index == 1 {
            target = TargetTranslation::new(subject, &code);
        }
        let f: Vec<_> = line.split('\t').collect();
        let frame: i8 = f[1].parse().unwrap();
        let offset: i32 = f[2].parse().unwrap();
        let end: i32 = f[3].parse().unwrap();
        let (_, stop, _) = target.get(frame, offset, end).unwrap();
        let context = if frame > 0 { frame - 1 } else { 2 - frame };
        assert_eq!(
            target.frames[context as usize].start,
            f[4].parse::<usize>().unwrap()
        );
        assert_eq!(stop, f[5].parse::<usize>().unwrap());
    }
    let expected: Vec<GappedHsp> = fs::read_to_string(format!("{root}gapped_hsps.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            GappedHsp {
                frame: f[9].parse().unwrap(),
                score: f[3].parse().unwrap(),
                q_start: f[6].parse().unwrap(),
                q_end: f[7].parse().unwrap(),
                q_gapped_start: f[8].parse().unwrap(),
                s_start: f[10].parse().unwrap(),
                s_end: f[11].parse().unwrap(),
                s_gapped_start: f[12].parse().unwrap(),
            }
        })
        .collect();
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&[query], subject, 1, profile).unwrap();
    let gapped = stages
        .gapped
        .iter()
        .map(|(_, hsp, _)| *hsp)
        .collect::<Vec<_>>();
    assert_eq!(gapped, expected);
    let actual =
        full_translation_traceback_blosum62(query, subject, 1, &gapped, 11, 1, 64).unwrap();
    let mut pass = 0;
    let mut oracle = Vec::new();
    for line in fs::read_to_string(format!("{root}traceback_events.tsv"))
        .unwrap()
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split('\t').collect();
        if f[0] == "TRACEBACK_INPUT" {
            pass += 1;
        }
        if pass == 2 && f[0] == "TRACEBACK_OUT_HSP" {
            oracle.push((
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i8>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            ));
        }
    }
    assert_eq!(pass, 2);
    let actual: Vec<_> = actual
        .iter()
        .map(|(hsp, retried)| {
            assert!(*retried);
            (
                hsp.score,
                hsp.frame,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
            )
        })
        .collect();
    assert_eq!(actual, oracle);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1154-1228:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:432:
// Only a fence causes full-subject retry.
#[test]
fn nonzero_base_partial_traceback_without_fence_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_hsp_no_fence_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let code = GeneticCode::try_from_id(1).unwrap();
    let mut target = TargetTranslation::new(subject, &code);
    for line in fs::read_to_string(format!("{root}translation_windows.tsv"))
        .unwrap()
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split('\t').collect();
        let frame: i8 = f[1].parse().unwrap();
        let offset: i32 = f[2].parse().unwrap();
        let end: i32 = f[3].parse().unwrap();
        let (_, stop, base) = target.get(frame, offset, end).unwrap();
        let context = if frame > 0 { frame - 1 } else { 2 - frame };
        assert_eq!(
            target.frames[context as usize].start,
            f[4].parse::<usize>().unwrap()
        );
        assert_eq!(stop, f[5].parse::<usize>().unwrap());
        assert!(base > 0);
    }
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 12,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&[query], subject, 1, profile).unwrap();
    let gapped = stages
        .gapped
        .iter()
        .map(|(_, hsp, _)| *hsp)
        .collect::<Vec<_>>();
    let oracle_gapped: Vec<_> = fs::read_to_string(format!("{root}gapped_hsps.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse::<i32>().unwrap(),
                f[9].parse::<i8>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[10].parse::<i32>().unwrap(),
                f[11].parse::<i32>().unwrap(),
            )
        })
        .collect();
    let actual_gapped: Vec<_> = gapped
        .iter()
        .map(|h| (h.score, h.frame, h.q_start, h.q_end, h.s_start, h.s_end))
        .collect();
    assert_eq!(actual_gapped, oracle_gapped);
    let actual =
        full_translation_traceback_blosum62(query, subject, 1, &gapped, 11, 1, 12).unwrap();
    let oracle: Vec<_> = fs::read_to_string(format!("{root}traceback_events.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .filter(|line| line.starts_with("TRACEBACK_OUT_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i8>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            )
        })
        .collect();
    let actual: Vec<_> = actual
        .iter()
        .map(|(h, retried)| {
            assert!(!retried);
            (h.score, h.frame, h.q_start, h.q_end, h.s_start, h.s_end)
        })
        .collect();
    assert_eq!(actual, oracle);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:804-844:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:635-693:
#[test]
fn six_frame_merged_traceback_endpoint_deletion_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/six_frame_merged_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let expected_gapped: Vec<GappedHsp> = fs::read_to_string(format!("{root}gapped_hsps.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            GappedHsp {
                frame: f[9].parse().unwrap(),
                score: f[3].parse().unwrap(),
                q_start: f[6].parse().unwrap(),
                q_end: f[7].parse().unwrap(),
                q_gapped_start: f[8].parse().unwrap(),
                s_start: f[10].parse().unwrap(),
                s_end: f[11].parse().unwrap(),
                s_gapped_start: f[12].parse().unwrap(),
            }
        })
        .collect();
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&[query], subject, 1, profile).unwrap();
    let gapped = stages
        .gapped
        .iter()
        .map(|(_, hsp, _)| *hsp)
        .collect::<Vec<_>>();
    assert_eq!(gapped, expected_gapped);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-589,635-666:
    // The extra -1 HSP grows from raw 38 to raw 236 before endpoint purge.
    let removed = gapped
        .iter()
        .find(|hsp| hsp.frame == -1 && hsp.score == 38)
        .unwrap();
    let code = GeneticCode::try_from_id(1).unwrap();
    let mut target = TargetTranslation::new(subject, &code);
    let (subject_sequence, _, _) = target.get(-1, -1, 0).unwrap();
    let encoded = encode_protein_query_frame_with_seg(query, None);
    let query_sequence = &encoded.aa_seq[1..encoded.aa_seq.len() - 1];
    let mut fence = false;
    let mut scratch = GapAlignScratch::new();
    let before_purge = blast_gapped_alignment_with_traceback_with_scratch(
        query_sequence,
        subject_sequence,
        removed.q_gapped_start as usize,
        removed.s_gapped_start as usize,
        ScoringMatrix::Blosum62,
        None,
        11,
        1,
        64,
        &mut scratch,
        Some(&mut fence),
    )
    .unwrap();
    assert!(!fence);
    let endpoint = fs::read_to_string(format!("{root}endpoint_events.tsv")).unwrap();
    let mut in_traceback_purge = false;
    let mut expected_removed = None;
    for line in endpoint.lines().skip(1) {
        if line == "ENDPOINT_INPUT\t0\t18" {
            in_traceback_purge = true;
        }
        if in_traceback_purge && line.starts_with("ENDPOINT_IN_HSP\t6\t") {
            let f: Vec<_> = line.split('\t').collect();
            expected_removed = Some((
                f[2].parse::<i32>().unwrap(),
                // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88: Int4 query_start/query_stop/subject_start/subject_stop.
                f[4].parse::<i32>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
            ));
            break;
        }
    }
    assert_eq!(
        Some((
            before_purge.score,
            before_purge.query_start,
            before_purge.query_stop,
            before_purge.subject_start,
            before_purge.subject_stop
        )),
        expected_removed
    );
    let actual =
        full_translation_traceback_blosum62(query, subject, 1, &gapped, 11, 1, 64).unwrap();
    let mut pass = 0;
    let mut expected = Vec::new();
    for line in fs::read_to_string(format!("{root}traceback_events.tsv"))
        .unwrap()
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split('\t').collect();
        if f[0] == "TRACEBACK_INPUT" {
            pass += 1;
        }
        if pass == 2 && f[0] == "TRACEBACK_OUT_HSP" {
            expected.push((
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i8>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            ));
        }
    }
    assert_eq!(pass, 2);
    assert_eq!(gapped.len(), 18);
    assert_eq!(expected.len(), 17);
    let actual: Vec<_> = actual
        .iter()
        .map(|(hsp, retried)| {
            assert!(*retried);
            (
                hsp.score,
                hsp.frame,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
            )
        })
        .collect();
    assert_eq!(actual, expected);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3924-3927,4057-4091:
// BLAST_GetGappedScore saves each qualifying HSP in frame/chunk order.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:539-586,840-850:
// purge, merge, and append follow each chunk's gapped extension.
#[test]
fn alternate_matrix_word2_gapped_and_append_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 16,
        window: 60,
        word_xdrop: &[21, 21, 0],
        word_cutoff: &[0, 0, i32::MAX],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum45,
        word_size: 2,
        gap_open: 14,
        gap_extend: 2,
        gap_xdrop: 53,
        gapped_cutoff: &[0, 0, i32::MAX],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&refs, subject, 1, profile).unwrap();
    let gapped = stages.gapped.clone();
    let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("GAPPED_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            assert_eq!(f[5], "0");
            (
                f[4].parse::<usize>().unwrap(),
                GappedHsp {
                    score: f[3].parse().unwrap(),
                    frame: f[9].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[10].parse().unwrap(),
                    s_end: f[11].parse().unwrap(),
                    s_gapped_start: f[12].parse().unwrap(),
                },
                0usize,
            )
        })
        .collect();
    if gapped != expected {
        let first = gapped
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(gapped.len().min(expected.len()));
        panic!(
            "first BLOSUM45/word2 gapped difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
            gapped.get(first), expected.get(first), gapped.len(), expected.len()
        );
    }

    let snapshots = stages.appended.clone();
    let append_trace = fs::read_to_string(format!("{root}append_events.tsv")).unwrap();
    let expected_snapshots: Vec<Vec<_>> = (0..6)
        .map(|call| {
            append_trace
                .lines()
                .filter(|line| line.starts_with(&format!("APPEND_OUT_HSP\t{call}\t")))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (
                        f[3].parse::<usize>().unwrap(),
                        f[4].parse::<i32>().unwrap(),
                        f[5].parse::<i8>().unwrap(),
                        f[6].parse::<i32>().unwrap(),
                        f[7].parse::<i32>().unwrap(),
                        f[8].parse::<i32>().unwrap(),
                        f[9].parse::<i32>().unwrap(),
                    )
                })
                .collect()
        })
        .collect();
    let actual_snapshots: Vec<Vec<_>> = snapshots
        .iter()
        .map(|snapshot| {
            snapshot
                .iter()
                .map(|(context, hsp)| {
                    (
                        *context,
                        hsp.score,
                        hsp.frame,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                    )
                })
                .collect()
        })
        .collect();
    assert_eq!(actual_snapshots, expected_snapshots);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:401-405,583-605:
// The pinned BLOSUM45/word-2 trace contains positive traceback-side
// containment, including the full-translation retry.
#[test]
fn alternate_matrix_traceback_positive_containment_matches_ncbi() {
    use std::collections::HashMap;
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let append_trace = fs::read_to_string(format!("{root}append_events.tsv")).unwrap();
    let mut appended = Vec::new();
    let mut merged_by_key = HashMap::new();
    for line in append_trace.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f[0] == "MERGE_OUT_HSP" {
            let context: usize = f[3].parse().unwrap();
            let hsp = GappedHsp {
                score: f[4].parse().unwrap(),
                frame: f[5].parse().unwrap(),
                q_start: f[6].parse().unwrap(),
                q_end: f[7].parse().unwrap(),
                q_gapped_start: f[8].parse().unwrap(),
                s_start: f[9].parse().unwrap(),
                s_end: f[10].parse().unwrap(),
                s_gapped_start: f[11].parse().unwrap(),
            };
            merged_by_key.insert(
                (
                    context,
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                ),
                hsp,
            );
        } else if f[0] == "APPEND_OUT_HSP" && f[1] == "5" {
            appended.push((
                f[3].parse::<usize>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[5].parse::<i8>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
                f[9].parse::<i32>().unwrap(),
            ));
        }
    }
    assert_eq!(appended.len(), 22);
    let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
    // the same natural query/subject inputs pass through WordFinder, gapped extension, chunk
    // purge, merge and append (BLOSUM45/word-2) before the positive traceback interval-tree
    // containment call.
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 16,
        window: 60,
        word_xdrop: &[21, 21, 0],
        word_cutoff: &[0, 0, i32::MAX],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum45,
        word_size: 2,
        gap_open: 14,
        gap_extend: 2,
        gap_xdrop: 53,
        gapped_cutoff: &[0, 0, i32::MAX],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let integrated = preliminary_hsps_of_letters(&refs, subject, 1, profile).unwrap();
    let complete_path = &integrated;
    let actual_appended: Vec<_> = complete_path
        .iter()
        .map(|(context, hsp)| {
            (
                *context,
                hsp.score,
                hsp.frame,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
            )
        })
        .collect();
    assert_eq!(actual_appended, appended);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:572-586,840-850;
    // c++/src/algo/blast/core/blast_hits.c:2809-2864:
    // append keeps each merged HSP's gapped starts. Check the exact
    // traceback input state, which APPEND_OUT_HSP does not print.
    for (context, hsp) in complete_path {
        let key = (
            *context,
            hsp.score,
            hsp.frame,
            hsp.q_start,
            hsp.q_end,
            hsp.s_start,
            hsp.s_end,
        );
        assert_eq!(*hsp, *merged_by_key.get(&key).unwrap());
    }
    let trace = fs::read_to_string(format!("{root}traceback_events.tsv")).unwrap();
    let mut expected: Vec<Vec<(i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); queries.len()];
    let mut expected_containment = vec![Vec::new(); queries.len()];
    let mut context = 0usize;
    let mut successful_pass = false;
    let mut in_traceback = false;
    let mut positive_containment = 0;
    for line in trace.lines() {
        let f: Vec<_> = line.split('\t').collect();
        match f[0] {
            "TRACEBACK_INPUT" => in_traceback = true,
            "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
            "CONTAINS" if in_traceback => {
                let result = f[1] == "1";
                positive_containment += usize::from(result);
                expected_containment[context].push((
                    result,
                    f[2].parse::<i8>().unwrap(),
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                ));
            }
            "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
            "TRACEBACK_OUT_HSP" if successful_pass => expected[context].push((
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
            )),
            _ => {}
        }
    }
    assert!(positive_containment > 0);
    for context in [1usize, 0] {
        let per_query: Vec<_> = complete_path
            .iter()
            .filter(|(index, _)| *index == context)
            .map(|(_, hsp)| *hsp)
            .collect();
        let mut containment_events = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events(
            &queries[context].1,
            subject,
            1,
            &per_query,
            ScoringMatrix::Blosum45,
            14,
            2,
            88,
            0.0,
            None,
            None,
            None,
            Some(&mut containment_events),
            None,
        )
        .unwrap();
        assert_eq!(
            containment_events, expected_containment[context],
            "NCBI containment predicate and call order, context {context}"
        );
        let actual: Vec<_> = actual
            .into_iter()
            .map(|(hsp, _)| {
                (
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
            .collect();
        assert_eq!(actual, expected[context], "context {context}");
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-552:
// GetGappedScore receives the chunk-local translated subject before
// Blast_HSPListAdjustOffsets and Blast_HSPListsMerge.
#[test]
fn long_subject_chunk_local_gapped_hsps_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let query = &read_fasta(&format!("{root}long_chunk_20260924/query.faa"))[0].1;
    let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
    let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
    let mut subject = Vec::with_capacity(15_000_962);
    for _ in 0..4_999_950 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    assert_eq!(subject.len(), 15_000_962);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590:
    // Retained long-fixture PARAM rows give cutoff 28.
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[28],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[28],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&[query], &subject, 1, profile).unwrap();
    let chunked = stages.gapped.clone();
    assert_eq!(
        chunked.iter().map(|row| row.2).collect::<Vec<_>>(),
        [0, 4_999_900]
    );
    let actual: Vec<_> = chunked.iter().map(|(_, hsp, _)| *hsp).collect();
    let trace = fs::read_to_string(format!("{root}long_chunk_20260924/gapped_events.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("GAPPED_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            GappedHsp {
                score: f[3].parse().unwrap(),
                q_start: f[6].parse().unwrap(),
                q_end: f[7].parse().unwrap(),
                q_gapped_start: f[8].parse().unwrap(),
                frame: f[9].parse().unwrap(),
                s_start: f[10].parse().unwrap(),
                s_end: f[11].parse().unwrap(),
                s_gapped_start: f[12].parse().unwrap(),
            }
        })
        .collect();
    assert_eq!(actual, expected);
    assert_eq!(actual.len(), 2);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:572-586:
    let mut second = actual[1];
    second.s_start += 4_999_900;
    second.s_end += 4_999_900;
    second.s_gapped_start += 4_999_900;
    let merged = merge_subject_chunk_hsps(
        vec![(0, actual[0])],
        vec![(0, second)],
        4_999_900,
        100,
        i32::MAX as usize,
        true,
    );
    let merged_expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("MERGE_OUT_HSP\t1\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse().unwrap(),
                GappedHsp {
                    score: f[4].parse().unwrap(),
                    frame: f[5].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[9].parse().unwrap(),
                    s_end: f[10].parse().unwrap(),
                    s_gapped_start: f[11].parse().unwrap(),
                },
            )
        })
        .collect();
    assert_eq!(merged, merged_expected);
    let snapshots = stages.appended.clone();
    assert_eq!(snapshots.iter().map(Vec::len).collect::<Vec<_>>(), [1; 6]);
    assert_eq!(snapshots.last().unwrap(), &merged_expected);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3924-3927:
// each long chunk uses its query-context cutoff during gapped scoring.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
// Both chunks merge, then the
// frame list is appended in score order for both valid query contexts.
#[test]
fn long_multi_query_chunk_hsps_and_append_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let fixture = format!("{root}long_multi_query_20260924/");
    let queries = read_fasta(&format!("{fixture}query.faa"));
    let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
    let plus1 = read_fasta(&format!("{root}run_20260923/subjects.fna"))
        .into_iter()
        .find(|(id, _)| id == "plus1")
        .unwrap()
        .1;
    let mut subject = Vec::with_capacity(15_000_962);
    for _ in 0..4_999_950 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(&plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    assert_eq!(subject.len(), 15_000_962);
    let cutoffs = fs::read_to_string(format!("{fixture}context_cutoffs.tsv")).unwrap();
    let gapped_cutoff_rows: Vec<_> = cutoffs
        .lines()
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            assert_eq!(f[0], "GAPPED_CONTEXT_CUTOFF");
            (
                f[1].parse::<usize>().unwrap(),
                f[2].parse::<usize>().unwrap(),
                f[4].parse::<i32>().unwrap(),
            )
        })
        .collect();
    assert_eq!(gapped_cutoff_rows.len(), 6);
    for call in 0..2 {
        for (context, value) in [28, 25, i32::MAX].into_iter().enumerate() {
            assert_eq!(
                gapped_cutoff_rows[call * 3 + context],
                (call, context, value)
            );
        }
    }
    let expected_cutoffs = [28, 25, i32::MAX];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:547-583:
    // the saved probe records 28, 25, INT4_MAX in all 12 chunk calls.
    let word_trace = fs::read_to_string(format!("{fixture}word_context_cutoffs.tsv")).unwrap();
    let word_rows: Vec<_> = word_trace
        .lines()
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            assert_eq!(f[0], "WORD_CONTEXT_CUTOFF");
            (
                f[1].parse::<usize>().unwrap(),
                f[2].parse::<usize>().unwrap(),
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[5].parse::<i32>().unwrap(),
            )
        })
        .collect();
    assert_eq!(word_rows.len(), 36);
    for call in 0..12 {
        for (context, (init, drop, cutoff)) in [(16, 16, 28), (16, 16, 25), (0, 0, i32::MAX)]
            .into_iter()
            .enumerate()
        {
            assert_eq!(
                word_rows[call * 3 + context],
                (call, context, init, drop, cutoff)
            );
        }
    }
    let word_cutoffs = [28, 25, i32::MAX];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200-234,547-614:
    // after per-context extension/cutoff, sorted chunk-local init HSPs
    // enter GetGappedScore; compare all rows before any gapped deletion.
    let wordfinder_trace = fs::read_to_string(format!("{fixture}frame_chunks.tsv")).unwrap();
    let expected_initial: Vec<_> = wordfinder_trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            InitHsp {
                frame: [1, 2, 3, -1, -2, -3][call / 2],
                chunk_offset: if call.is_multiple_of(2) { 0 } else { 4_999_900 },
                q_seed: f[3].parse().unwrap(),
                s_seed: f[4].parse().unwrap(),
                q_start: f[5].parse().unwrap(),
                s_start: f[6].parse().unwrap(),
                length: f[7].parse().unwrap(),
                score: f[8].parse().unwrap(),
            }
        })
        .collect();
    let trace = fs::read_to_string(format!("{fixture}gapped_events.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("GAPPED_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            (
                f[4].parse::<usize>().unwrap(),
                GappedHsp {
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    frame: f[9].parse().unwrap(),
                    s_start: f[10].parse().unwrap(),
                    s_end: f[11].parse().unwrap(),
                    s_gapped_start: f[12].parse().unwrap(),
                },
                if call == 0 { 0 } else { 4_999_900 },
            )
        })
        .collect();
    let expected_snapshots: Vec<Vec<_>> = (0..6)
        .map(|call| {
            trace
                .lines()
                .filter(|line| line.starts_with(&format!("APPEND_OUT_HSP\t{call}\t")))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (
                        f[3].parse::<usize>().unwrap(),
                        f[4].parse::<i32>().unwrap(),
                        f[5].parse::<i8>().unwrap(),
                        f[6].parse::<i32>().unwrap(),
                        f[7].parse::<i32>().unwrap(),
                        f[8].parse::<i32>().unwrap(),
                        f[9].parse::<i32>().unwrap(),
                    )
                })
                .collect()
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:572-586,840-850:
    // The final merged HSPs
    // retain exact gapped starts through append into traceback.
    let expected_merged: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("MERGE_OUT_HSP\t1\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse::<usize>().unwrap(),
                GappedHsp {
                    score: f[4].parse().unwrap(),
                    frame: f[5].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[9].parse().unwrap(),
                    s_end: f[10].parse().unwrap(),
                    s_gapped_start: f[11].parse().unwrap(),
                },
            )
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
    // each chunk's scansub -> WordFinder -> GetGappedScore -> purge ->
    // offset adjustment -> merge precedes the next chunk; append ends
    // each frame. Compare this one-call Rust path with every saved NCBI
    // function input and output, including all 394 candidate pairs.
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16, 16, 0],
        word_cutoff: &word_cutoffs,
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &expected_cutoffs,
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) = preliminary_hsps_traced(&refs, &subject, 1, profile, None).unwrap();
    assert_eq!(integrated, expected_merged);
    let frame_order = [1, 2, 3, -1, -2, -3];
    let candidate_trace = fs::read_to_string(format!("{fixture}candidate.stderr")).unwrap();
    let mut candidate_by_call = vec![Vec::new(); 12];
    for line in candidate_trace
        .lines()
        .filter(|line| line.starts_with("CAND\t"))
    {
        let f: Vec<_> = line.split('\t').collect();
        let call: usize = f[1].parse().unwrap();
        candidate_by_call[call].push(Seed {
            frame: frame_order[call / 2],
            chunk_offset: if call.is_multiple_of(2) { 0 } else { 4_999_900 },
            query_offset: f[3].parse().unwrap(),
            subject_offset: f[4].parse().unwrap(),
        });
    }
    assert_eq!(candidate_by_call.iter().map(Vec::len).sum::<usize>(), 394);
    // The per-event compares below filter these lists by (frame, chunk); the sizes show that
    // no expected row is left outside the events.
    assert_eq!(expected_initial.len(), 4);
    assert_eq!(expected.len(), 4);
    let mut tags = Vec::new();
    let mut call = 0usize;
    let mut append = 0usize;
    for event in &events {
        match event {
            PreliminaryEvent::Candidates(frame, offset, seeds) => {
                tags.push("C");
                assert_eq!(
                    (*frame, *offset),
                    (
                        frame_order[call / 2],
                        if call.is_multiple_of(2) { 0 } else { 4_999_900 }
                    )
                );
                assert_eq!(seeds, &candidate_by_call[call], "candidate call {call}");
                call += 1;
            }
            PreliminaryEvent::Initial(frame, offset, hsps) => {
                tags.push("I");
                let expected: Vec<_> = expected_initial
                    .iter()
                    .copied()
                    .filter(|hit| hit.frame == *frame && hit.chunk_offset as usize == *offset)
                    .collect();
                assert_eq!(*hsps, expected, "WordFinder call {}", call - 1);
            }
            PreliminaryEvent::Gapped(frame, offset, hsps) => {
                tags.push("G");
                let expected: Vec<_> = expected
                    .iter()
                    .filter(|(_, hsp, chunk_offset)| hsp.frame == *frame && chunk_offset == offset)
                    .map(|(context, hsp, _)| (*context, *hsp))
                    .collect();
                assert_eq!(*hsps, expected, "GetGappedScore call {}", call - 1);
            }
            PreliminaryEvent::Purged(frame, offset, hsps) => {
                tags.push("P");
                let expected: Vec<_> = expected
                    .iter()
                    .filter(|(_, hsp, chunk_offset)| hsp.frame == *frame && chunk_offset == offset)
                    .map(|(context, hsp, _)| (*context, *hsp))
                    .collect();
                assert_eq!(*hsps, expected, "natural endpoint purge call {}", call - 1);
            }
            PreliminaryEvent::Merged(frame, offset, hsps) => {
                tags.push("M");
                let merge_call = if *offset == 0 { 0 } else { 1 };
                assert_eq!(*frame, 1);
                let expected: Vec<_> = trace
                    .lines()
                    .filter(|line| line.starts_with(&format!("MERGE_OUT_HSP\t{merge_call}\t")))
                    .map(|line| {
                        let f: Vec<_> = line.split('\t').collect();
                        (
                            f[3].parse::<usize>().unwrap(),
                            GappedHsp {
                                score: f[4].parse().unwrap(),
                                frame: f[5].parse().unwrap(),
                                q_start: f[6].parse().unwrap(),
                                q_end: f[7].parse().unwrap(),
                                q_gapped_start: f[8].parse().unwrap(),
                                s_start: f[9].parse().unwrap(),
                                s_end: f[10].parse().unwrap(),
                                s_gapped_start: f[11].parse().unwrap(),
                            },
                        )
                    })
                    .collect();
                assert_eq!(*hsps, expected, "merge call {merge_call}");
            }
            PreliminaryEvent::Appended(frame, hsps) => {
                tags.push("A");
                assert_eq!(*frame, frame_order[append]);
                let actual: Vec<_> = hsps
                    .iter()
                    .map(|(context, hsp)| {
                        (
                            *context,
                            hsp.score,
                            hsp.frame,
                            hsp.q_start,
                            hsp.q_end,
                            hsp.s_start,
                            hsp.s_end,
                        )
                    })
                    .collect();
                assert_eq!(actual, expected_snapshots[append], "append call {append}");
                append += 1;
            }
        }
    }
    let mut expected_tags = Vec::new();
    for call in 0..12 {
        expected_tags.extend(["C", "I"]);
        if call < 2 {
            expected_tags.extend(["G", "P", "M"]);
        }
        if call % 2 == 1 {
            expected_tags.push("A");
        }
    }
    assert_eq!(tags, expected_tags, "NCBI function execution order");
    let count = |pick: fn(&PreliminaryEvent) -> Option<usize>| -> usize {
        events.iter().filter_map(pick).sum()
    };
    assert_eq!(
        count(|event| match event {
            PreliminaryEvent::Initial(_, _, hsps) => Some(hsps.len()),
            _ => None,
        }),
        4
    );
    assert_eq!(
        count(|event| match event {
            PreliminaryEvent::Gapped(_, _, hsps) => Some(hsps.len()),
            _ => None,
        }),
        4
    );
    assert_eq!((call, append), (12, 6));
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:259-312,401-449,583-721:
    // Process each
    // query-indexed appended list; a fence retries on full translation.
    let mut expected_final = vec![Vec::new(); refs.len()];
    let mut context = 0usize;
    let mut successful_pass = false;
    for line in trace.lines() {
        let f: Vec<_> = line.split('\t').collect();
        match f[0] {
            "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
            "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
            "TRACEBACK_OUT_HSP" if successful_pass => expected_final[context].push((
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i8>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            )),
            _ => {}
        }
    }
    for context in [1usize, 0] {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:259-312:
        // each query-indexed HSP list enters traceback from the completed
        // frame append stream; the saved trace visits context 1 then 0.
        let per_query: Vec<_> = integrated
            .iter()
            .filter(|(index, _)| *index == context)
            .map(|(_, hsp)| *hsp)
            .collect();
        let actual =
            full_translation_traceback_blosum62(refs[context], &subject, 1, &per_query, 11, 1, 64)
                .unwrap();
        let actual: Vec<_> = actual
            .into_iter()
            .map(|(hsp, retry)| {
                assert!(retry);
                (
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
            .collect();
        assert_eq!(actual, expected_final[context], "query context {context}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:283-310,478-586:
// the masked middle chunk returns SUBJECT_SPLIT_NO_RANGE; the later
// positive-frame chunk yields one local HSP, then gets offset-adjusted.
#[test]
fn no_range_middle_chunk_gapped_and_merged_hsp_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let query = &read_fasta(&format!("{root}no_range_middle_20260924/query.faa"))[0].1;
    let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
    let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
    let mut subject = Vec::with_capacity(30_001_262);
    for _ in 0..10_000_050 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    subject[14_997_000..30_000_000].make_ascii_lowercase();
    let trace =
        fs::read_to_string(format!("{root}no_range_middle_20260924/gapped_events.tsv")).unwrap();
    let local: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("GAPPED_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[4].parse().unwrap(),
                GappedHsp {
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    frame: f[9].parse().unwrap(),
                    s_start: f[10].parse().unwrap(),
                    s_end: f[11].parse().unwrap(),
                    s_gapped_start: f[12].parse().unwrap(),
                },
                9_999_800,
            )
        })
        .collect();
    let merged: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("MERGE_OUT_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse().unwrap(),
                GappedHsp {
                    score: f[4].parse().unwrap(),
                    frame: f[5].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[9].parse().unwrap(),
                    s_end: f[10].parse().unwrap(),
                    s_gapped_start: f[11].parse().unwrap(),
                },
            )
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590:
    // Retained long-fixture PARAM rows give cutoff 30.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
    // SUBJECT_SPLIT_NO_RANGE skips WordFinder/GetGappedScore in the
    // masked middle chunk; append still occurs once for every frame.
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[30],
        mask_lowercase: true,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[30],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) =
        preliminary_hsps_traced(&[query], &subject, 1, profile, None).unwrap();
    assert_eq!(integrated, merged);
    // The gapped lists and append snapshots of the events (offsets ride in the events).
    let chunked: Vec<_> = events
        .iter()
        .filter_map(|event| match event {
            PreliminaryEvent::Gapped(_, offset, hsps) => Some((*offset, hsps)),
            _ => None,
        })
        .flat_map(|(offset, hsps)| {
            hsps.iter()
                .map(move |(context, hsp)| (*context, *hsp, offset))
        })
        .collect();
    assert_eq!(
        chunked.iter().map(|row| row.2).collect::<Vec<_>>(),
        [9_999_800]
    );
    assert_eq!(chunked, local);
    let snapshots: Vec<_> = events
        .iter()
        .filter_map(|event| match event {
            PreliminaryEvent::Appended(_, hsps) => Some(hsps.clone()),
            _ => None,
        })
        .collect();
    assert_eq!(snapshots.iter().map(Vec::len).collect::<Vec<_>>(), [1; 6]);
    assert_eq!(snapshots.last().unwrap(), &merged);
    assert_eq!(
        events
            .iter()
            .filter(|event| matches!(event, PreliminaryEvent::Appended(..)))
            .count(),
        6
    );
    assert_eq!(
        events
            .iter()
            .filter(|event| matches!(event, PreliminaryEvent::Gapped(..)))
            .count(),
        1
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-625;
// c++/src/algo/blast/core/blast_engine.c:484-525;
// c++/src/algo/blast/core/blast_gapalign.c:2410-2442:
// the hard-SEG working query is shared by WordFinder and GetGappedScore;
// the unmasked suffix still produces a real HSP and traceback.
#[test]
fn hard_seg_query_bytes_and_integrated_hsp_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/seg_hard_query_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let seg = SegParams::default();
    let encoded = encode_protein_query_frame_with_seg(query, Some(&seg));
    assert!(!encoded.seg_masks.is_empty());
    assert!(encoded.aa_seq[1..41].iter().all(|&b| b == 21));
    let hex: String = encoded.aa_seq[1..1 + query.len()]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
    let word: Vec<_> = state
        .lines()
        .filter(|line| line.starts_with("WORD_QUERY_BYTES\t"))
        .collect();
    let gap: Vec<_> = state
        .lines()
        .filter(|line| line.starts_with("GAPPED_QUERY_BYTES\t"))
        .collect();
    assert_eq!((word.len(), gap.len()), (6, 1));
    for (call, row) in word.iter().enumerate() {
        assert_eq!(*row, format!("WORD_QUERY_BYTES\t{call}\t160\t{hex}"));
    }
    assert_eq!(gap[0], format!("GAPPED_QUERY_BYTES\t0\t160\t{hex}"));
    let profile = PreliminaryProfile {
        seg: Some(&seg),
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) =
        preliminary_hsps_traced(&[query], subject, 1, profile, None).unwrap();
    let candidate_trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
    let mut candidates = vec![Vec::new(); 6];
    for line in candidate_trace
        .lines()
        .filter(|line| line.starts_with("CAND\t"))
    {
        let f: Vec<_> = line.split('\t').collect();
        let call: usize = f[1].parse().unwrap();
        candidates[call].push(Seed {
            frame: [1, 2, 3, -1, -2, -3][call],
            chunk_offset: 0,
            query_offset: f[3].parse().unwrap(),
            subject_offset: f[4].parse().unwrap(),
        });
    }
    assert_eq!(candidates.iter().map(Vec::len).sum::<usize>(), 185);
    let gapped_trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let fields: Vec<_> = gapped_trace
        .lines()
        .find(|line| line.starts_with("GAPPED_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    let ncbi_hsp = GappedHsp {
        score: fields[3].parse().unwrap(),
        frame: fields[9].parse().unwrap(),
        q_start: fields[6].parse().unwrap(),
        q_end: fields[7].parse().unwrap(),
        q_gapped_start: fields[8].parse().unwrap(),
        s_start: fields[10].parse().unwrap(),
        s_end: fields[11].parse().unwrap(),
        s_gapped_start: fields[12].parse().unwrap(),
    };
    assert_eq!(integrated, vec![(0, ncbi_hsp)]);
    let mut call = 0usize;
    let mut append = 0usize;
    for event in &events {
        match event {
            PreliminaryEvent::Candidates(frame, offset, seeds) => {
                assert_eq!(
                    (*frame, *offset, seeds),
                    ([1, 2, 3, -1, -2, -3][call], 0, &candidates[call])
                );
                call += 1;
            }
            PreliminaryEvent::Initial(frame, _, hsps) => {
                assert_eq!(*frame, [1, 2, 3, -1, -2, -3][call - 1]);
                assert_eq!(hsps.len(), usize::from(call == 1));
                if call == 1 {
                    assert_eq!(
                        (
                            hsps[0].q_seed,
                            hsps[0].s_seed,
                            hsps[0].q_start,
                            hsps[0].s_start,
                            hsps[0].length,
                            hsps[0].score
                        ),
                        (43, 3, 40, 0, 120, 656)
                    );
                }
            }
            PreliminaryEvent::Gapped(_, _, hsps) | PreliminaryEvent::Purged(_, _, hsps) => {
                assert_eq!(*hsps, integrated)
            }
            PreliminaryEvent::Merged(_, _, hsps) => assert_eq!(*hsps, integrated),
            PreliminaryEvent::Appended(frame, hsps) => {
                assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                assert_eq!(*hsps, integrated);
                append += 1;
            }
        }
    }
    assert_eq!((call, append), (6, 6));
    let traceback =
        full_translation_traceback_blosum62(query, subject, 1, &[integrated[0].1], 11, 1, 64)
            .unwrap();
    assert_eq!(traceback, vec![(ncbi_hsp, true)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:380-391,583-596:
// query_blk->sequence (masked) drives the traceback score, while
// query_blk->sequence_nomask drives identities before Blast_HSPTest.
// The retained NCBI call has a real HSP spanning masked residues 60..78.
#[test]
fn internal_seg_crossing_traceback_uses_masked_score_and_unmasked_identity() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/seg_cross_traceback_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let seg = SegParams::default();
    let encoded = encode_protein_query_frame_with_seg(query, Some(&seg));
    assert_eq!(query.len(), 138);
    assert!(encoded.aa_seq[61..79].iter().all(|&b| b == 21));
    let unmasked = encoded.aa_seq_nomask.as_ref().unwrap();
    assert!(unmasked[61..79].iter().all(|&b| b == 10));
    let masked_hex: String = encoded.aa_seq[1..139]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let unmasked_hex: String = unmasked[1..139]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
    let word: Vec<_> = state
        .lines()
        .filter(|line| line.starts_with("WORD_QUERY_BYTES\t"))
        .collect();
    let gap: Vec<_> = state
        .lines()
        .filter(|line| line.starts_with("GAPPED_QUERY_BYTES\t"))
        .collect();
    let identity: Vec<_> = state
        .lines()
        .filter(|line| line.starts_with("IDENTITY_STATE\t"))
        .collect();
    assert_eq!((word.len(), gap.len(), identity.len()), (6, 1, 1));
    for (call, row) in word.iter().enumerate() {
        assert_eq!(*row, format!("WORD_QUERY_BYTES\t{call}\t138\t{masked_hex}"));
    }
    assert_eq!(gap[0], format!("GAPPED_QUERY_BYTES\t0\t138\t{masked_hex}"));
    assert_eq!(
        identity[0],
        format!("IDENTITY_STATE\t0\t638\t138\t138\t0\t138\t{unmasked_hex}")
    );
    let profile = PreliminaryProfile {
        seg: Some(&seg),
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) =
        preliminary_hsps_traced(&[query], subject, 1, profile, None).unwrap();
    let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let f: Vec<_> = trace
        .lines()
        .find(|line| line.starts_with("GAPPED_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    let expected = GappedHsp {
        score: f[3].parse().unwrap(),
        frame: f[9].parse().unwrap(),
        q_start: f[6].parse().unwrap(),
        q_end: f[7].parse().unwrap(),
        q_gapped_start: f[8].parse().unwrap(),
        s_start: f[10].parse().unwrap(),
        s_end: f[11].parse().unwrap(),
        s_gapped_start: f[12].parse().unwrap(),
    };
    assert_eq!(expected.score, 638);
    assert_eq!((expected.q_start, expected.q_end), (0, 138));
    assert_eq!(integrated, vec![(0, expected)]);
    assert_eq!(
        events
            .iter()
            .filter(|event| matches!(event, PreliminaryEvent::Appended(..)))
            .count(),
        6
    );
    let mut identities = Vec::new();
    let actual = full_translation_traceback_with_matrix_and_events(
        query,
        subject,
        1,
        &[integrated[0].1],
        ScoringMatrix::Blosum62,
        11,
        1,
        64,
        0.0,
        Some(&seg),
        None,
        None,
        None,
        Some(&mut identities),
    )
    .unwrap();
    assert_eq!(actual, vec![(expected, true)]);
    assert_eq!(identities, [(138, 138)]);
    let unmasked_result = full_translation_traceback_with_matrix(
        query,
        subject,
        1,
        &[integrated[0].1],
        ScoringMatrix::Blosum62,
        11,
        1,
        64,
        0.0,
    )
    .unwrap();
    assert_eq!(unmasked_result[0].0.score, 746);
    assert_ne!(unmasked_result[0].0.score, expected.score);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-638;
// c++/src/algo/blast/core/blast_engine.c:484-525;
// c++/src/algo/blast/core/blast_traceback.c:380-391,583-596:
// Soft masking excludes the same seed intervals but keeps K in every
// WordFinder, GetGappedScore and traceback query input.
#[test]
fn internal_soft_seg_crossing_matches_ncbi_query_hsp_and_traceback() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/seg_soft_traceback_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let seg = SegParams::default();
    let unmasked = encode_protein_query_frame_with_seg(query, None);
    let hex: String = unmasked.aa_seq[1..139]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
    let word: Vec<_> = state
        .lines()
        .filter(|line| line.starts_with("WORD_QUERY_BYTES\t"))
        .collect();
    assert_eq!(word.len(), 6);
    for (call, row) in word.iter().enumerate() {
        assert_eq!(*row, format!("WORD_QUERY_BYTES\t{call}\t138\t{hex}"));
    }
    assert!(state
        .lines()
        .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t138\t{hex}")));
    assert!(state
        .lines()
        .any(|row| row == format!("IDENTITY_STATE\t0\t746\t138\t138\t0\t138\t{hex}")));
    let profile = PreliminaryProfile {
        seg: Some(&seg),
        soft_masking: true,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: false,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) =
        preliminary_hsps_traced(&[query], subject, 1, profile, None).unwrap();
    let candidate_trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
    let mut candidates = vec![Vec::new(); 6];
    for line in candidate_trace
        .lines()
        .filter(|line| line.starts_with("CAND\t"))
    {
        let f: Vec<_> = line.split('\t').collect();
        let call: usize = f[1].parse().unwrap();
        candidates[call].push(Seed {
            frame: [1, 2, 3, -1, -2, -3][call],
            chunk_offset: 0,
            query_offset: f[3].parse().unwrap(),
            subject_offset: f[4].parse().unwrap(),
        });
    }
    assert_eq!(candidates.iter().map(Vec::len).sum::<usize>(), 183);
    let mut call = 0;
    let mut append = 0;
    for event in &events {
        match event {
            PreliminaryEvent::Candidates(frame, offset, seeds) => {
                assert_eq!(
                    (*frame, *offset, seeds),
                    ([1, 2, 3, -1, -2, -3][call], 0, &candidates[call])
                );
                call += 1;
            }
            PreliminaryEvent::Initial(_, _, hsps) => {
                assert_eq!(hsps.len(), usize::from(call == 1));
                if call == 1 {
                    assert_eq!(
                        (
                            hsps[0].q_seed,
                            hsps[0].s_seed,
                            hsps[0].q_start,
                            hsps[0].s_start,
                            hsps[0].length,
                            hsps[0].score
                        ),
                        (3, 3, 0, 0, 138, 746)
                    );
                }
            }
            PreliminaryEvent::Appended(frame, hsps) => {
                assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                assert_eq!(*hsps, integrated);
                append += 1;
            }
            _ => {}
        }
    }
    assert_eq!((call, append), (6, 6));
    let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let f: Vec<_> = trace
        .lines()
        .find(|line| line.starts_with("GAPPED_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    let expected = GappedHsp {
        score: f[3].parse().unwrap(),
        frame: f[9].parse().unwrap(),
        q_start: f[6].parse().unwrap(),
        q_end: f[7].parse().unwrap(),
        q_gapped_start: f[8].parse().unwrap(),
        s_start: f[10].parse().unwrap(),
        s_end: f[11].parse().unwrap(),
        s_gapped_start: f[12].parse().unwrap(),
    };
    assert_eq!(expected.score, 746);
    assert_eq!(integrated, vec![(0, expected)]);
    let mut identities = Vec::new();
    let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
        query,
        subject,
        1,
        &[expected],
        ScoringMatrix::Blosum62,
        11,
        1,
        64,
        0.0,
        Some(&seg),
        true,
        false,
        None,
        None,
        None,
        Some(&mut identities),
    )
    .unwrap();
    assert_eq!(actual, vec![(expected, true)]);
    assert_eq!(identities, [(138, 138)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:486-502;
// c++/src/algo/blast/core/blast_filter.c:1241-1255;
// c++/src/algo/blast/core/blast_setup.c:614-638:
// -lcase_masking includes query lowercase spans in the filter locations.
// The hard-mask working query has X in positions 40..60; identities use
// the original sequence after the traceback score is computed.
#[test]
fn lowercase_query_mask_matches_ncbi_candidates_hsps_and_traceback() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/lcase_query_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let masked = encode_tblastn_lookup_query(query, None, true);
    assert_eq!(query.len(), 120);
    assert!(masked.aa_seq[41..61].iter().all(|&residue| residue == 21));
    let nomask = masked.aa_seq_nomask.as_ref().unwrap();
    let masked_hex: String = masked.aa_seq[1..121]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let nomask_hex: String = nomask[1..121]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
    for call in 0..6 {
        assert!(state
            .lines()
            .any(|row| row == format!("WORD_QUERY_BYTES\t{call}\t120\t{masked_hex}")));
    }
    assert!(state
        .lines()
        .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t120\t{masked_hex}")));
    assert!(state
        .lines()
        .any(|row| row == format!("IDENTITY_STATE\t0\t528\t120\t120\t0\t120\t{nomask_hex}")));
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: true,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) =
        preliminary_hsps_traced(&[query], subject, 1, profile, None).unwrap();
    let candidate_trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
    let mut candidates = vec![Vec::new(); 6];
    for line in candidate_trace
        .lines()
        .filter(|line| line.starts_with("CAND\t"))
    {
        let f: Vec<_> = line.split('\t').collect();
        let call: usize = f[1].parse().unwrap();
        candidates[call].push(Seed {
            frame: [1, 2, 3, -1, -2, -3][call],
            chunk_offset: 0,
            query_offset: f[3].parse().unwrap(),
            subject_offset: f[4].parse().unwrap(),
        });
    }
    assert_eq!(candidates.iter().map(Vec::len).sum::<usize>(), 154);
    let initial_trace = fs::read_to_string(format!("{root}wordfinder.stderr")).unwrap();
    let initial: Vec<_> = initial_trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse::<u32>().unwrap(),
                f[4].parse::<u32>().unwrap(),
                f[5].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
            )
        })
        .collect();
    assert_eq!(initial, [(63, 63, 60, 60, 60, 325), (3, 3, 0, 0, 40, 223)]);
    let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let f: Vec<_> = trace
        .lines()
        .find(|line| line.starts_with("GAPPED_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    let expected = GappedHsp {
        score: f[3].parse().unwrap(),
        frame: f[9].parse().unwrap(),
        q_start: f[6].parse().unwrap(),
        q_end: f[7].parse().unwrap(),
        q_gapped_start: f[8].parse().unwrap(),
        s_start: f[10].parse().unwrap(),
        s_end: f[11].parse().unwrap(),
        s_gapped_start: f[12].parse().unwrap(),
    };
    assert_eq!(expected.score, 528);
    assert_eq!(integrated, vec![(0, expected)]);
    let mut call = 0;
    let mut append = 0;
    for event in &events {
        match event {
            PreliminaryEvent::Candidates(frame, offset, seeds) => {
                assert_eq!(
                    (*frame, *offset, seeds),
                    ([1, 2, 3, -1, -2, -3][call], 0, &candidates[call])
                );
                call += 1;
            }
            PreliminaryEvent::Initial(_, _, hsps) => {
                assert_eq!(
                    hsps.iter()
                        .map(|h| (h.q_seed, h.s_seed, h.q_start, h.s_start, h.length, h.score))
                        .collect::<Vec<_>>(),
                    if call == 1 {
                        initial.clone()
                    } else {
                        Vec::new()
                    }
                );
            }
            PreliminaryEvent::Appended(frame, hsps) => {
                assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                assert_eq!(*hsps, integrated);
                append += 1;
            }
            _ => {}
        }
    }
    assert_eq!((call, append), (6, 6));
    let mut identities = Vec::new();
    let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
        query,
        subject,
        1,
        &[expected],
        ScoringMatrix::Blosum62,
        11,
        1,
        64,
        0.0,
        None,
        false,
        true,
        None,
        None,
        None,
        Some(&mut identities),
    )
    .unwrap();
    assert_eq!(actual, vec![(expected, true)]);
    assert_eq!(identities, [(120, 120)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-638;
// c++/src/algo/blast/core/blast_filter.c:1241-1255:
// mask_at_hash leaves the lowercase query unmasked during extension and
// traceback, while its complement still restricts lookup positions.
#[test]
fn lowercase_query_soft_mask_matches_ncbi_lookup_and_traceback() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/lcase_soft_query_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-638 keeps this query unmasked when
    // mask_at_hash is true; all six WordFinder calls see the same bytes.
    let unmasked = encode_protein_query_frame_with_seg(query, None);
    let hex: String = unmasked.aa_seq[1..121]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
    for call in 0..6 {
        assert!(state
            .lines()
            .any(|row| row == format!("WORD_QUERY_BYTES\t{call}\t120\t{hex}")));
    }
    assert!(state
        .lines()
        .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t120\t{hex}")));
    assert!(state
        .lines()
        .any(|row| row == format!("IDENTITY_STATE\t0\t656\t120\t120\t0\t120\t{hex}")));
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: true,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: true,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) =
        preliminary_hsps_traced(&[query], subject, 1, profile, None).unwrap();
    let reference = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let f: Vec<_> = reference
        .lines()
        .find(|line| line.starts_with("GAPPED_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    let expected = GappedHsp {
        score: f[3].parse().unwrap(),
        frame: f[9].parse().unwrap(),
        q_start: f[6].parse().unwrap(),
        q_end: f[7].parse().unwrap(),
        q_gapped_start: f[8].parse().unwrap(),
        s_start: f[10].parse().unwrap(),
        s_end: f[11].parse().unwrap(),
        s_gapped_start: f[12].parse().unwrap(),
    };
    assert_eq!(expected.score, 656);
    assert_eq!(integrated, vec![(0, expected)]);
    let hard = fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/lcase_query_20260924/candidate.stderr"
    ))
    .unwrap();
    assert_eq!(
        fs::read_to_string(format!("{root}candidate.stderr")).unwrap(),
        hard
    );
    let mut saved_candidates = vec![Vec::new(); 6];
    for line in hard.lines().filter(|line| line.starts_with("CAND\t")) {
        let f: Vec<_> = line.split('\t').collect();
        let call: usize = f[1].parse().unwrap();
        saved_candidates[call].push(Seed {
            frame: [1, 2, 3, -1, -2, -3][call],
            chunk_offset: 0,
            query_offset: f[3].parse().unwrap(),
            subject_offset: f[4].parse().unwrap(),
        });
    }
    let mut candidate_call = 0;
    let mut initial_call = 0;
    let mut append_call = 0;
    for event in &events {
        match event {
            PreliminaryEvent::Candidates(frame, offset, seeds) => {
                assert_eq!(
                    (*frame, *offset, seeds),
                    (
                        [1, 2, 3, -1, -2, -3][candidate_call],
                        0,
                        &saved_candidates[candidate_call]
                    )
                );
                candidate_call += 1;
            }
            PreliminaryEvent::Initial(_, _, hsps) => {
                assert_eq!(hsps.len(), usize::from(initial_call == 0));
                if initial_call == 0 {
                    assert_eq!(
                        (
                            hsps[0].q_seed,
                            hsps[0].s_seed,
                            hsps[0].q_start,
                            hsps[0].s_start,
                            hsps[0].length,
                            hsps[0].score
                        ),
                        (3, 3, 0, 0, 120, 656)
                    );
                }
                initial_call += 1;
            }
            PreliminaryEvent::Appended(frame, hsps) => {
                assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append_call]);
                assert_eq!(*hsps, integrated);
                append_call += 1;
            }
            _ => {}
        }
    }
    assert_eq!((candidate_call, initial_call, append_call), (6, 6, 6));
    let mut identities = Vec::new();
    let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
        query,
        subject,
        1,
        &[expected],
        ScoringMatrix::Blosum62,
        11,
        1,
        64,
        0.0,
        None,
        true,
        true,
        None,
        None,
        None,
        Some(&mut identities),
    )
    .unwrap();
    assert_eq!(actual, vec![(expected, true)]);
    assert_eq!(identities, [(120, 120)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1241-1256;
// c++/src/algo/blast/core/blast_setup.c:614-638:
// The SEG prefix and overlapping lowercase span combine into [0,45).
#[test]
fn overlapping_seg_and_lowercase_query_mask_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/seg_lcase_overlap_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let seg = SegParams::default();
    let masked = encode_tblastn_lookup_query(query, Some(&seg), true);
    assert_eq!(masked.seg_masks, [(0, 45)]);
    assert!(masked.aa_seq[1..46].iter().all(|&residue| residue == 21));
    let hex: String = masked.aa_seq[1..161]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
    for call in 0..6 {
        assert!(state
            .lines()
            .any(|row| row == format!("WORD_QUERY_BYTES\t{call}\t160\t{hex}")));
    }
    assert!(state
        .lines()
        .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t160\t{hex}")));
    let profile = PreliminaryProfile {
        seg: Some(&seg),
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: true,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) =
        preliminary_hsps_traced(&[query], subject, 1, profile, None).unwrap();
    let candidate_trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
    let mut saved_candidates = vec![Vec::new(); 6];
    for line in candidate_trace
        .lines()
        .filter(|line| line.starts_with("CAND\t"))
    {
        let f: Vec<_> = line.split('\t').collect();
        let call: usize = f[1].parse().unwrap();
        saved_candidates[call].push(Seed {
            frame: [1, 2, 3, -1, -2, -3][call],
            chunk_offset: 0,
            query_offset: f[3].parse().unwrap(),
            subject_offset: f[4].parse().unwrap(),
        });
    }
    assert_eq!(saved_candidates.iter().map(Vec::len).sum::<usize>(), 168);
    let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let f: Vec<_> = trace
        .lines()
        .find(|line| line.starts_with("GAPPED_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    let expected = GappedHsp {
        score: f[3].parse().unwrap(),
        frame: f[9].parse().unwrap(),
        q_start: f[6].parse().unwrap(),
        q_end: f[7].parse().unwrap(),
        q_gapped_start: f[8].parse().unwrap(),
        s_start: f[10].parse().unwrap(),
        s_end: f[11].parse().unwrap(),
        s_gapped_start: f[12].parse().unwrap(),
    };
    assert_eq!(expected.score, 623);
    assert_eq!(integrated, vec![(0, expected)]);
    let mut call = 0;
    let mut append = 0;
    for event in &events {
        match event {
            PreliminaryEvent::Candidates(frame, offset, seeds) => {
                assert_eq!(
                    (*frame, *offset, seeds),
                    ([1, 2, 3, -1, -2, -3][call], 0, &saved_candidates[call])
                );
                call += 1;
            }
            PreliminaryEvent::Initial(_, _, hsps) => {
                assert_eq!(hsps.len(), usize::from(call == 1));
                if call == 1 {
                    assert_eq!(
                        (
                            hsps[0].q_seed,
                            hsps[0].s_seed,
                            hsps[0].q_start,
                            hsps[0].s_start,
                            hsps[0].length,
                            hsps[0].score
                        ),
                        (48, 8, 45, 5, 115, 623)
                    );
                }
            }
            PreliminaryEvent::Appended(frame, hsps) => {
                assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                assert_eq!(*hsps, integrated);
                append += 1;
            }
            _ => {}
        }
    }
    assert_eq!((call, append), (6, 6));
    let mut identities = Vec::new();
    let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
        query,
        subject,
        1,
        &[expected],
        ScoringMatrix::Blosum62,
        11,
        1,
        64,
        0.0,
        Some(&seg),
        false,
        true,
        None,
        None,
        None,
        Some(&mut identities),
    )
    .unwrap();
    assert_eq!(actual, vec![(expected, true)]);
    assert_eq!(identities, [(115, 115)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1241-1256;
// c++/src/algo/blast/core/blast_setup.c:614-638:
// The overlapping SEG/lowercase intervals still restrict lookup, while
// WordFinder, GetGappedScore and traceback use the unmasked query.
#[test]
fn overlapping_seg_lowercase_soft_mask_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/seg_lcase_overlap_soft_20260924/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let seg = SegParams::default();
    let unmasked = encode_protein_query_frame_with_seg(query, None);
    let hex: String = unmasked.aa_seq[1..161]
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
    for call in 0..6 {
        assert!(state
            .lines()
            .any(|row| row == format!("WORD_QUERY_BYTES\t{call}\t160\t{hex}")));
    }
    assert!(state
        .lines()
        .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t160\t{hex}")));
    assert!(state
        .lines()
        .any(|row| row == format!("IDENTITY_STATE\t0\t656\t120\t120\t40\t160\t{hex}")));
    let profile = PreliminaryProfile {
        seg: Some(&seg),
        soft_masking: true,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[0],
        mask_lowercase: true,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[0],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let (integrated, events) =
        preliminary_hsps_traced(&[query], subject, 1, profile, None).unwrap();
    let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
    let f: Vec<_> = trace
        .lines()
        .find(|line| line.starts_with("GAPPED_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    let expected = GappedHsp {
        score: f[3].parse().unwrap(),
        frame: f[9].parse().unwrap(),
        q_start: f[6].parse().unwrap(),
        q_end: f[7].parse().unwrap(),
        q_gapped_start: f[8].parse().unwrap(),
        s_start: f[10].parse().unwrap(),
        s_end: f[11].parse().unwrap(),
        s_gapped_start: f[12].parse().unwrap(),
    };
    assert_eq!(expected.score, 656);
    assert_eq!(integrated, vec![(0, expected)]);
    let candidates = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
    let hard_candidates = fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/seg_lcase_overlap_20260924/candidate.stderr"
    ))
    .unwrap();
    assert_eq!(candidates, hard_candidates);
    let mut saved_candidates = vec![Vec::new(); 6];
    for line in candidates.lines().filter(|line| line.starts_with("CAND\t")) {
        let f: Vec<_> = line.split('\t').collect();
        let call: usize = f[1].parse().unwrap();
        saved_candidates[call].push(Seed {
            frame: [1, 2, 3, -1, -2, -3][call],
            chunk_offset: 0,
            query_offset: f[3].parse().unwrap(),
            subject_offset: f[4].parse().unwrap(),
        });
    }
    assert_eq!(saved_candidates.iter().map(Vec::len).sum::<usize>(), 168);
    let mut call = 0;
    let mut append = 0;
    for event in &events {
        match event {
            PreliminaryEvent::Candidates(frame, offset, seeds) => {
                assert_eq!(
                    (*frame, *offset, seeds),
                    ([1, 2, 3, -1, -2, -3][call], 0, &saved_candidates[call])
                );
                call += 1;
            }
            PreliminaryEvent::Initial(_, _, hsps) => {
                assert_eq!(hsps.len(), usize::from(call == 1));
                if call == 1 {
                    assert_eq!(
                        (
                            hsps[0].q_seed,
                            hsps[0].s_seed,
                            hsps[0].q_start,
                            hsps[0].s_start,
                            hsps[0].length,
                            hsps[0].score
                        ),
                        (48, 8, 40, 0, 120, 656)
                    );
                }
            }
            PreliminaryEvent::Appended(frame, hsps) => {
                assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                assert_eq!(*hsps, integrated);
                append += 1;
            }
            _ => {}
        }
    }
    assert_eq!((call, append), (6, 6));
    let mut identities = Vec::new();
    let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
        query,
        subject,
        1,
        &[expected],
        ScoringMatrix::Blosum62,
        11,
        1,
        64,
        0.0,
        Some(&seg),
        true,
        true,
        None,
        None,
        None,
        Some(&mut identities),
    )
    .unwrap();
    assert_eq!(actual, vec![(expected, true)]);
    assert_eq!(identities, [(120, 120)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
// skip the fully masked middle chunk; retain the first chunk's list,
// then merge the later chunk and append that frame before the next frame.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-3035:
// Blast_HSPListsMerge retains both nonoverlapping equal-score HSPs.
#[test]
fn skipped_middle_chunk_preserves_both_hsps_and_all_append_rows() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let case = format!("{root}no_range_two_hits_20260924/");
    let query = &read_fasta(&format!("{case}query.faa"))[0].1;
    let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
    let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
    let mut subject = Vec::with_capacity(30_001_262);
    for _ in 0..4_998_450 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(plus1);
    subject.push(b'A');
    for _ in 0..(10_000_050 - 4_998_450 - 121) {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    assert_eq!(subject.len(), 30_001_262);
    subject[14_997_000..30_000_000].make_ascii_lowercase();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590:
    // Retained long-fixture PARAM rows give cutoff 30.
    let profile = PreliminaryProfile {
        seg: None,
        soft_masking: false,
        threshold: 13,
        window: 40,
        word_xdrop: &[16],
        word_cutoff: &[30],
        mask_lowercase: true,
        matrix: ScoringMatrix::Blosum62,
        word_size: 3,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &[30],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    };
    let stages = preliminary_stages(&[query], &subject, 1, profile).unwrap();
    let initial = stages.initial.clone();
    let word_trace = fs::read_to_string(format!("{case}frame_chunks.tsv")).unwrap();
    let expected_init: Vec<_> = word_trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            assert!(call < 2);
            super::super::search_init::InitHsp {
                frame: 1,
                chunk_offset: if call == 0 { 0 } else { 9_999_800 },
                q_seed: f[3].parse().unwrap(),
                s_seed: f[4].parse().unwrap(),
                q_start: f[5].parse().unwrap(),
                s_start: f[6].parse().unwrap(),
                length: f[7].parse().unwrap(),
                score: f[8].parse().unwrap(),
            }
        })
        .collect();
    assert_eq!(initial, expected_init);
    let chunked = stages.gapped.clone();
    let trace = fs::read_to_string(format!("{case}gapped_events.tsv")).unwrap();
    let expected_local: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("GAPPED_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            assert!(call < 2);
            (
                f[4].parse().unwrap(),
                GappedHsp {
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    frame: f[9].parse().unwrap(),
                    s_start: f[10].parse().unwrap(),
                    s_end: f[11].parse().unwrap(),
                    s_gapped_start: f[12].parse().unwrap(),
                },
                if call == 0 { 0 } else { 9_999_800 },
            )
        })
        .collect();
    assert_eq!(chunked, expected_local);
    let snapshots = stages.appended.clone();
    let expected_append: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("APPEND_OUT_HSP\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].parse::<usize>().unwrap(),
                f[3].parse::<usize>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[5].parse::<i8>().unwrap(),
                f[6].parse::<i32>().unwrap(),
                f[7].parse::<i32>().unwrap(),
                f[8].parse::<i32>().unwrap(),
                f[9].parse::<i32>().unwrap(),
            )
        })
        .collect();
    let actual_append: Vec<_> = snapshots
        .iter()
        .enumerate()
        .flat_map(|(call, hsps)| {
            hsps.iter().map(move |(context, hsp)| {
                (
                    call,
                    *context,
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
        })
        .collect();
    assert_eq!(actual_append, expected_append);
    assert_eq!(snapshots.iter().map(Vec::len).collect::<Vec<_>>(), [2; 6]);
}

// The HSPs of an ungapped unit (`BLAST_GetUngappedHSPList`, cited at `ungapped_chunk_hsps`):
// each initial HSP in the context of its word's query offset (a word at a context's first
// offset is that context's), the query offsets made the context's, the word's offsets as the
// gapped start, the subject offsets the chunk's, and the list sorted by score.
#[test]
fn ungapped_chunk_hsps_are_context_local_and_sorted_by_score() {
    // Two queries of 10 and 13 letters: contexts at 0 and 11 of the concatenated query.
    let contexts = [0, 11];
    let hit = |q_seed: u32, q_start: i32, s_seed: u32, s_start: i32, length, score| InitHsp {
        frame: 2,
        chunk_offset: 0,
        q_seed,
        s_seed,
        q_start,
        s_start,
        length,
        score,
    };
    let hsps = ungapped_chunk_hsps(
        &contexts,
        2,
        &[
            hit(15, 13, 102, 100, 6, 30),
            hit(3, 1, 42, 40, 5, 50),
            hit(11, 11, 7, 7, 3, 30),
        ],
    );
    let hsp = |score, q_start, q_end, q_gapped_start, s_start, s_end, s_gapped_start| GappedHsp {
        frame: 2,
        score,
        q_start,
        q_end,
        q_gapped_start,
        s_start,
        s_end,
        s_gapped_start,
    };
    assert_eq!(
        hsps,
        [
            (0, hsp(50, 1, 6, 3, 40, 45, 42)),
            // Equal scores: the smaller subject offset first (`ScoreCompareHSPs`).
            (1, hsp(30, 0, 3, 0, 7, 10, 7)),
            (1, hsp(30, 2, 8, 4, 100, 106, 102)),
        ]
    );
}
