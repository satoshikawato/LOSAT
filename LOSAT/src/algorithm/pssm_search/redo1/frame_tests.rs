//! Tests of the redo's windows of a translated subject (rpstblastn, item PS12b): NCBI's list of
//! windows (`s_WindowsFromTranslatedAligns` with a position-based matrix), the frame copy of
//! `s_SequenceGetProteinRange` in a buffer of the translation buffer's length plus 2, and the
//! composition range of a window that ends before it begins (`-strand plus|minus`).

use super::*;
use crate::algorithm::pssm_search::gapalign::PssmGapScratch;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::stats::KarlinParams;
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::utils::matrix::{aa_char_to_ncbistdaa, ncbistdaa};

fn hsp(frame: i32, s: (i32, i32), score: i32) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: 0,
        query_end: s.1 - s.0,
        subject_start: s.0,
        subject_end: s.1,
        gapped_query_start: 0,
        gapped_subject_start: s.0,
        raw_score: score,
        query_context: 0,
        query_frame: 0,
        subject_frame: frame,
        query_length: 0,
        q_idx: 0,
        s_idx: 0,
    }
}

fn windows(
    hsps: &[BlastpPreliminaryHsp],
    subject_length: i32,
) -> Vec<((i32, i32, i32), Vec<usize>)> {
    frame_windows(hsps, 20, subject_length)
        .into_iter()
        .map(|(w, order)| ((w.begin, w.end, w.frame), order))
        .collect()
}

/// A redo run of the query `consensus` (only the near-identical test reads it).
fn run_of(consensus: &[u8]) -> RedoRun<'_> {
    RedoRun::new(
        consensus,
        consensus,
        KarlinParams {
            lambda: 0.267 / SCALING_FACTOR,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        },
        0.041f64.glibc_ln(),
        (11, 1),
        (25.0, 0),
        QueryPssm {
            rows: vec![[0; crate::core::blast_util::BLASTAA_SIZE]; consensus.len()],
        },
        PssmGapScratch::new(),
    )
}

fn encode(text: &str) -> Vec<u8> {
    text.bytes().map(aa_char_to_ncbistdaa).collect()
}

// redo_alignment.c:644-765 (cited at `windows_from_translated_aligns`): one window per alignment,
// `[max(0, matchStart - 200), min(translated length, matchEnd + 200))` in its frame, with the
// translated length of `GET_NUCL_LENGTH` of the buffer's length; sorted by query, frame, begin
// and end (glibc's stable sort), a window joined into the next one of its frame that it reaches
// (the next one's alignments first), each window's alignments by score, then the windows by
// begin, end and frame.
#[test]
fn the_windows_are_ncbis() {
    // A buffer of 121 bytes (both strands of 60 nucleotides): frames +1 and -1 hold 20 letters,
    // +2 19. Frame +1's two alignments share a window, the better one first; the window of
    // frame +2 ends first, then frame -1 sorts before frame +1.
    let hsps = [
        hsp(1, (2, 10), 50),
        hsp(-1, (0, 8), 40),
        hsp(1, (12, 18), 60),
        hsp(2, (5, 15), 45),
    ];
    assert_eq!(
        windows(&hsps, 121),
        [
            ((0, 19, 2), vec![3]),
            ((0, 20, -1), vec![1]),
            ((0, 20, 1), vec![2, 0])
        ]
    );
    // 2,001 bytes (1,000 nucleotides, frame +1: 333 letters): windows that overlap are joined;
    // equal scores sort by the alignments' starts.
    let hsps = [hsp(1, (250, 260), 30), hsp(1, (10, 20), 30)];
    assert_eq!(windows(&hsps, 2001), [((0, 333, 1), vec![1, 0])]);
    // 4,001 bytes (frame +1: 666 letters): windows apart stay apart.
    let hsps = [hsp(1, (450, 460), 30), hsp(1, (10, 20), 30)];
    assert_eq!(
        windows(&hsps, 4001),
        [((0, 220, 1), vec![1]), ((250, 660, 1), vec![0])]
    );
}

// blast_kappa.c:1597-1643 (cited at `frame_range`): the frame's letters as NCBI finds them in
// the translation buffer, in a zeroed buffer of the buffer's length plus 2, the byte before the
// window zeroed.
#[test]
fn the_frame_copy_is_ncbis() {
    let consensus = encode("MKVLWAALLV");
    let run = run_of(&consensus);
    let pattern = |length: usize| -> Vec<u8> { (0..length).map(|i| (i % 20 + 1) as u8).collect() };
    let mut biased = false;
    // Both strands (121 bytes): frame -1 starts after three frames and their sentinels (21 +
    // 20 + 20 bytes).
    let subject = pattern(121);
    let frame = align_in(-1, (6, 12));
    let buffer = frame_range(
        &run,
        &subject,
        window(5, 20, -1),
        &frame,
        false,
        &mut biased,
    );
    assert_eq!(buffer.len(), 123);
    let mut expected = vec![0u8; 123];
    expected[1..21].copy_from_slice(&subject[61..81]);
    expected[5] = 0;
    assert_eq!(buffer, expected);
    // -strand plus (61 bytes: three frames of 20, 19, 19 letters and their sentinels): NCBI's
    // arithmetic reads 30 nucleotides, frame -3 as the 9 bytes from byte 52, the buffer's end.
    let subject = pattern(61);
    let frame = align_in(-3, (0, 4));
    let buffer = frame_range(&run, &subject, window(0, 9, -3), &frame, false, &mut biased);
    assert_eq!(buffer.len(), 63);
    assert_eq!(&buffer[1..10], &subject[52..61]);
    assert!(buffer[10..].iter().all(|&b| b == 0) && buffer[0] == 0);
    // SEG masks the frame's copy (not the buffer) where the subject may be biased: a frame of
    // one letter is masked, a frame of 20 different letters is not (no longer biased).
    let mut subject = encode(&format!("{}*{}", "ACDEFGHIKLMNPQRSTVWY", "A".repeat(19)));
    subject.resize(121, ncbistdaa::A);
    let mut biased = true;
    let buffer = frame_range(
        &run,
        &subject,
        window(0, 19, 2),
        &align_in(2, (0, 5)),
        false,
        &mut biased,
    );
    assert!(biased);
    assert!(buffer[1..20].iter().all(|&b| b == ncbistdaa::X));
    assert_eq!(buffer[20], 0);
    let buffer = frame_range(
        &run,
        &subject,
        window(0, 20, 1),
        &align_in(1, (0, 5)),
        false,
        &mut biased,
    );
    assert!(!biased);
    assert_eq!(&buffer[1..21], &subject[..20]);
}

fn window(begin: i32, end: i32, frame: i32) -> SubjectWindow {
    SubjectWindow { begin, end, frame }
}

fn align_in(frame: i32, s: (i32, i32)) -> Align {
    Align::incoming::<FrameRedo>(&hsp(frame, s, 30))
}

// redo_alignment.c:921-936 and composition_adjustment.c:1204-1248 (`Blast_GetCompositionRange`,
// cited at `blast_get_composition_range`): the letters between the stops around the alignment,
// 20 letters inside each stop unless that cuts into the alignment, within the window's letters;
// where the window ends before it begins, the range ends at the alignment's end.
#[test]
fn the_window_composition_is_ncbis() {
    let mut data = vec![ncbistdaa::A; 302];
    // A window of 150 letters from 0, stops at its letters 50 and 145: the alignment [100, 120)
    // reads [71, 125).
    data[1 + 50] = ncbistdaa::STOP;
    data[1 + 145] = ncbistdaa::STOP;
    let composition = window_composition(&data, window(0, 150, 1), &align_in(1, (100, 120)));
    assert_eq!(composition.num_true_amino_acids, 54);
    // A stop 10 letters after the alignment: the range ends at the alignment's end.
    let composition = window_composition(&data, window(0, 150, 1), &align_in(1, (100, 135)));
    assert_eq!(composition.num_true_amino_acids, 135 - 71);
    // A window [50, 40) (`-strand plus`: the alignment [250, 260) past the frame's 40 letters):
    // its letters read from byte 1 + 50, a stop at letter 100 of the window: [121, 210).
    data[1 + 50 + 100] = ncbistdaa::STOP;
    let composition = window_composition(&data, window(50, 40, 1), &align_in(1, (250, 260)));
    assert_eq!(composition.num_true_amino_acids, 89);
}

// redo_alignment.c:1063-1094 (cited at `preliminary_test_near_identical` and
// `old_test_near_identical`): with the cutoff 0 (rpstblastn's query whose frame +1 is invalid,
// `RedoRun::without_valid_context`) the window must hold one ungapped alignment as long as the
// profile or 50, whatever its score; with a cutoff above 0 the score per position decides.
#[test]
fn the_old_near_identical_test_is_ncbis() {
    let long = align_in(-1, (10, 70));
    assert!(old_test_near_identical(200, 1, &long));
    assert!(!old_test_near_identical(200, 2, &long));
    assert!(!old_test_near_identical(200, 1, &align_in(-1, (10, 40))));
    assert!(old_test_near_identical(25, 1, &align_in(-1, (10, 40))));
    let mut gapped = hsp(-1, (10, 70), 30);
    gapped.query_end += 1;
    assert!(!old_test_near_identical(
        200,
        1,
        &Align::incoming::<FrameRedo>(&gapped)
    ));
    // The same alignment with the cutoff of a valid frame +1 (1.74 bits per position): 960 / 60.
    assert!(!preliminary_test_near_identical(200, &long, 144.548446));
}
