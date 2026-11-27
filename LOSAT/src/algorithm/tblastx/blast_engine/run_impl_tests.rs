use super::*;

fn make_chunk_hit(
    ctx_idx: usize,
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    score: i32,
) -> UngappedHit {
    UngappedHit {
        q_idx: 0,
        s_idx: 0,
        ctx_idx: u32::try_from(ctx_idx).unwrap(),
        s_f_idx: 0,
        q_frame: 1,
        s_frame: 1,
        q_aa_start: u32::try_from(q_start).unwrap(),
        q_aa_end: u32::try_from(q_end).unwrap(),
        s_aa_start: u32::try_from(s_start).unwrap(),
        s_aa_end: u32::try_from(s_end).unwrap(),
        q_seed_off: u32::try_from(q_start).unwrap(),
        s_seed_off: u32::try_from(s_start).unwrap(),
        q_orig_len: 900,
        s_orig_len: 900,
        raw_score: score,
        e_value: f64::INFINITY,
        num_ident: 0,
        hsp_list_order: 0,
        ordering_method: 0,
        linked_set: false,
        start_of_chain: false,
        link_id: 0,
        chain_next_link_id: None,
        hsp_link_num: 0,
        num: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-310
#[test]
fn subject_split_state_uses_ncbi_max_len_and_overlap() {
    let mut split = SubjectSplitState::new(MAX_DBSEQ_LEN + 250);

    assert_eq!(
        split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
        SubjectChunkStatus::Ok(SubjectChunk {
            offset: 0,
            length: MAX_DBSEQ_LEN,
            overlap: 0,
            ncbi_length: MAX_DBSEQ_LEN as i32,
            overlap_past_list: false,
        })
    );
    assert_eq!(
        split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
        SubjectChunkStatus::Ok(SubjectChunk {
            offset: MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP,
            length: DBSEQ_CHUNK_OVERLAP + 250,
            overlap: DBSEQ_CHUNK_OVERLAP,
            ncbi_length: (DBSEQ_CHUNK_OVERLAP + 250) as i32,
            overlap_past_list: false,
        })
    );
    assert_eq!(
        split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
        SubjectChunkStatus::Done
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-2995
#[test]
fn merge_tblastx_subject_chunk_hits_merges_overlap_strip_same_context() {
    let mut combined = vec![make_chunk_hit(0, 0, 120, 0, 120, 100)];
    let incoming = vec![make_chunk_hit(0, 100, 220, 100, 220, 120)];

    merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP)
        .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED);

    assert_eq!(combined.len(), 1);
    assert_eq!(combined[0].q_aa_start, 0);
    assert_eq!(combined[0].q_aa_end, 220);
    assert_eq!(combined[0].s_aa_start, 0);
    assert_eq!(combined[0].s_aa_end, 220);
    assert_eq!(combined[0].q_seed_off, 100);
    assert_eq!(combined[0].s_seed_off, 100);
    assert_eq!(combined[0].raw_score, 201);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2968-2970
#[test]
fn merge_tblastx_subject_chunk_hits_keeps_different_context_hits() {
    let mut combined = vec![make_chunk_hit(0, 0, 120, 0, 120, 100)];
    let incoming = vec![make_chunk_hit(1, 100, 220, 100, 220, 120)];

    merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP)
        .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED);

    assert_eq!(combined.len(), 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-310
#[test]
fn subject_split_state_keeps_short_subjects_as_single_chunk() {
    let mut split = SubjectSplitState::new(MAX_DBSEQ_LEN - 1);

    assert_eq!(
        split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
        SubjectChunkStatus::Ok(SubjectChunk {
            offset: 0,
            length: MAX_DBSEQ_LEN - 1,
            overlap: 0,
            ncbi_length: (MAX_DBSEQ_LEN - 1) as i32,
            overlap_past_list: false,
        })
    );
    assert_eq!(
        split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
        SubjectChunkStatus::Done
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 578-584
// A hard-masked frame takes the chunks of its ranges in order, each with its merge overlap;
// a chunk that runs onto the end sentinel of a frame one residue shorter than its strand's
// first frame is cut at the frame's end.
#[test]
fn a_hard_masked_frame_lists_its_chunks_cut_at_the_frame() {
    let ranges = [(0, 300), (390, 1200)];
    let frame = |frame: i8, aa_len: usize| QueryFrame {
        frame,
        aa_seq: Vec::new(),
        aa_seq_nomask: None,
        aa_len,
        orig_len: 1200,
        seg_masks: Vec::new(),
    };
    let frames = [frame(1, 400), frame(2, 399), frame(3, 399), frame(-1, 400)];
    let chunks = |state: &mut SubjectSplitState| {
        std::iter::from_fn(
            || match state.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP) {
                SubjectChunkStatus::Ok(chunk) => Some((chunk.offset, chunk.length, chunk.overlap)),
                SubjectChunkStatus::Done => None,
            },
        )
        .collect::<Vec<_>>()
    };
    let calls = x_db_mask_calls(&ranges, true, &frames[1], &frames);
    assert_eq!(calls[1].0.length, 270, "the frame's ranges end one past it");
    let mut state = SubjectSplitState::listed(399, &calls);
    assert_eq!(chunks(&mut state), [(0, 100, 100), (130, 269, 100)]);
    let calls = x_db_mask_calls(&ranges, true, &frames[3], &frames);
    let mut state = SubjectSplitState::listed(400, &calls);
    assert_eq!(chunks(&mut state), [(0, 270, 100), (300, 100, 100)]);
    // With soft masks the frame's chunk is the whole frame, scanned in its ranges.
    let calls = x_db_mask_calls(&ranges, false, &frames[0], &frames);
    assert_eq!(calls.len(), 1);
    assert_eq!((calls[0].0.offset, calls[0].0.length), (0, 400));
    assert_eq!(calls[0].1, [(0, 100), (130, 400)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-2995
#[test]
fn merge_tblastx_subject_chunk_hits_keeps_previous_hits_before_overlap_strip() {
    let mut combined = vec![make_chunk_hit(0, 0, 90, 0, 90, 100)];
    let incoming = vec![make_chunk_hit(0, 100, 180, 100, 180, 120)];

    merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP)
        .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED);

    assert_eq!(combined.len(), 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-2995
#[test]
fn merge_tblastx_subject_chunk_hits_keeps_current_hits_after_overlap_strip() {
    let mut combined = vec![make_chunk_hit(0, 0, 140, 0, 140, 100)];
    let incoming = vec![make_chunk_hit(0, 200, 260, 200, 260, 120)];

    merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP)
        .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED);

    assert_eq!(combined.len(), 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1487-1531
#[test]
fn merge_tblastx_subject_chunk_hits_merges_hsp_fully_inside_overlap_strip() {
    let mut combined = vec![make_chunk_hit(0, 0, 190, 0, 190, 190)];
    let incoming = vec![make_chunk_hit(0, 120, 180, 120, 180, 60)];

    merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP)
        .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED);

    assert_eq!(combined.len(), 1);
    assert_eq!(combined[0].q_aa_start, 0);
    assert_eq!(combined[0].q_aa_end, 190);
    assert_eq!(combined[0].s_aa_start, 0);
    assert_eq!(combined[0].s_aa_end, 190);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1487-1531
#[test]
fn merge_tblastx_subject_chunk_hits_merges_negative_subject_frame_overlap() {
    let mut combined = vec![make_chunk_hit(0, 0, 120, 0, 120, 100)];
    combined[0].s_frame = -1;
    let mut incoming_hit = make_chunk_hit(0, 100, 220, 100, 220, 120);
    incoming_hit.s_frame = -1;

    merge_tblastx_subject_chunk_hits(&mut combined, vec![incoming_hit], 100, DBSEQ_CHUNK_OVERLAP)
        .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED);

    assert_eq!(combined.len(), 1);
    assert_eq!(combined[0].s_frame, -1);
    assert_eq!(combined[0].q_aa_end, 220);
    assert_eq!(combined[0].s_aa_end, 220);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3038-3051
#[test]
fn adjust_tblastx_chunk_subject_offsets_adjusts_subject_coordinates() {
    let mut hits = vec![make_chunk_hit(0, 10, 40, 20, 50, 100)];

    adjust_tblastx_chunk_subject_offsets(&mut hits, 500);

    assert_eq!(hits[0].s_aa_start, 520);
    assert_eq!(hits[0].s_aa_end, 550);
    assert_eq!(hits[0].s_seed_off, 520);
    assert_eq!(hits[0].q_aa_start, 10);
    assert_eq!(hits[0].q_seed_off, 10);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:83-127
#[test]
fn scan_interior_ranges_keep_right_lookahead_without_duplicate_emission() {
    let wordsize = 3usize;
    let subject_len = 8usize;
    let base_ranges = [(0, subject_len as i32)];

    let left = clip_tblastx_seq_ranges_for_scan_interior(&base_ranges, 0, 4, wordsize, subject_len);
    let right =
        clip_tblastx_seq_ranges_for_scan_interior(&base_ranges, 4, 8, wordsize, subject_len);

    assert_eq!(left, vec![(0, 6)]);
    assert_eq!(right, vec![(4, 8)]);

    let left_emitted: Vec<i32> = (left[0].0..=left[0].1 - wordsize as i32).collect();
    let right_emitted: Vec<i32> = (right[0].0..=right[0].1 - wordsize as i32).collect();
    assert_eq!(left_emitted, vec![0, 1, 2, 3]);
    assert_eq!(right_emitted, vec![4, 5]);
}

/// `make_chunk_hit` with NCBI's `Int4` ends of an HSP of the one-hit word finder.
fn make_int4_chunk_hit(
    ctx_idx: usize,
    q_start: usize,
    q_end: i32,
    s_start: usize,
    s_end: i32,
    score: i32,
) -> UngappedHit {
    let mut hit = make_chunk_hit(ctx_idx, q_start, q_start, s_start, s_start, score);
    hit.set_int4_ends(q_end, s_end);
    hit
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3037-3051
// An end below 0 (an HSP of length -2 at the start of a chunk) moves as an `Int4`.
#[test]
fn chunk_offsets_move_an_end_below_zero() {
    let mut hits = vec![
        make_int4_chunk_hit(0, 1, -1, 0, -2, 7),
        make_chunk_hit(0, 10, 40, 20, 50, 30),
    ];
    adjust_tblastx_chunk_subject_offsets(&mut hits, 500);
    assert_eq!(hits[0].s_aa_start, 500);
    assert_eq!(hits[0].s_end_int4(), 498);
    assert_eq!(hits[0].s_seed_off, 500);
    assert_eq!(hits[0].q_end_int4(), -1);
    assert!(hits[0].is_int4_inverted());
    assert_eq!(
        (hits[1].s_aa_start, hits[1].s_aa_end, hits[1].s_seed_off),
        (520, 550, 520)
    );
    assert_eq!((hits[1].q_aa_start, hits[1].q_aa_end), (10, 40));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3037-3051
// An HSP of the chunk from residue 4,999,900 that starts at its offset -95, and a seed at -1 (masks
// stored out of order, M4b-T2), move as NCBI's `Int4`.
#[test]
fn chunk_offsets_move_a_start_and_a_seed_before_the_chunk() {
    let mut hit = make_chunk_hit(0, 0, 60, 0, 0, 100);
    hit.s_aa_start = (-95i32) as u32;
    hit.set_int4_ends(60, -35);
    hit.s_seed_off = (-1i32) as u32;
    let mut hits = vec![hit];
    adjust_tblastx_chunk_subject_offsets(&mut hits, 4_999_900);
    assert_eq!(
        (hits[0].s_aa_start, hits[0].s_end_int4(), hits[0].s_seed_off),
        (4_999_805, 4_999_865, 4_999_899)
    );
    assert_eq!((hits[0].q_aa_start, hits[0].q_end_int4()), (0, 60));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1513-1527
// A positive HSP of the earlier chunk absorbs an HSP of length -1 that starts inside it: its
// bounds stay, and the density (100 + 7) / (20 - 1) raises its score to trunc(107 / 19 * 20) = 112.
#[test]
fn chunk_merge_absorbs_an_hsp_of_negative_length() {
    let mut combined = vec![make_chunk_hit(0, 100, 120, 100, 120, 100)];
    let incoming = vec![make_int4_chunk_hit(0, 110, 109, 110, 109, 7)];
    assert!(merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, 100).is_ok());
    assert_eq!(combined.len(), 1);
    let hit = &combined[0];
    assert_eq!((hit.q_aa_start, hit.q_end_int4()), (100, 120));
    assert_eq!((hit.s_aa_start, hit.s_end_int4()), (100, 120));
    assert_eq!(hit.raw_score, (107.0f64 / 19.0 * 20.0) as i32);
    assert_eq!(hit.raw_score, 112);
    assert_eq!((hit.q_seed_off, hit.s_seed_off), (100, 100));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits_priv.h:68-69
// An HSP of length 0 of the earlier chunk contains the HSP that starts at its point: it takes
// that HSP's bounds and seed (a higher score), and its score becomes trunc((5 + 30) / (0 + 10) *
// 10) = 35. An HSP of negative length contains nothing.
#[test]
fn chunk_merge_lets_an_hsp_of_length_zero_absorb_and_one_of_negative_length_not() {
    let mut combined = vec![make_int4_chunk_hit(0, 110, 110, 120, 120, 5)];
    let incoming = vec![make_chunk_hit(0, 110, 120, 120, 130, 30)];
    assert!(merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, 100).is_ok());
    assert_eq!(combined.len(), 1);
    let hit = &combined[0];
    assert!(!hit.is_int4_inverted());
    assert_eq!((hit.q_aa_start, hit.q_end_int4()), (110, 120));
    assert_eq!((hit.s_aa_start, hit.s_end_int4()), (120, 130));
    assert_eq!(hit.raw_score, 35);
    assert_eq!((hit.q_seed_off, hit.s_seed_off), (110, 120));

    let mut combined = vec![make_int4_chunk_hit(0, 110, 109, 120, 119, 5)];
    let incoming = vec![make_chunk_hit(0, 110, 120, 120, 130, 30)];
    assert!(merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, 100).is_ok());
    assert_eq!(combined.len(), 2);
    assert!(combined.iter().any(UngappedHit::is_int4_inverted));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1513-1515,1517-1520
// Lengths that add up to 0 (an HSP of length 0 found by both chunks, or one of length 2 and
// one of length -2): NCBI's density is +inf and `(int)` of its product undefined, refused
// (rule 3). The HSPs of the later chunk are still put in the list, which the refused search
// never reports.
#[test]
fn chunk_merge_refuses_lengths_that_add_up_to_zero() {
    let mut combined = vec![make_int4_chunk_hit(0, 110, 110, 120, 120, 5)];
    let incoming = vec![make_int4_chunk_hit(0, 110, 110, 120, 120, 5)];
    assert!(merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, 100).is_err());
    assert_eq!(combined.len(), 2);

    let mut combined = vec![make_chunk_hit(0, 110, 112, 120, 122, 9)];
    let incoming = vec![make_int4_chunk_hit(0, 111, 109, 121, 119, 4)];
    assert!(merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, 100).is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2893-2899
// Two HSPs of the earlier chunk that tie on every key of the score sort (two contexts, the
// same coordinates and score) ahead of a lower one in the overlap strip: the swap moves the
// strip HSP to the front and the first of the two behind the second ([N1, N2, S] becomes
// [S, N2, N1]), and the stable sort after the append keeps N2 before N1.
#[test]
fn chunk_merge_keeps_ncbis_order_of_tied_hsps() {
    let mut combined = vec![
        make_chunk_hit(0, 10, 20, 10, 20, 30),
        make_chunk_hit(1, 10, 20, 10, 20, 30),
        make_chunk_hit(2, 190, 250, 190, 250, 10),
    ];
    let incoming = vec![make_chunk_hit(0, 300, 320, 310, 330, 5)];
    for (order, hit) in combined.iter_mut().enumerate() {
        hit.hsp_list_order = order as u32;
    }
    assert!(merge_tblastx_subject_chunk_hits(&mut combined, incoming, 200, 100).is_ok());
    let contexts: Vec<u32> = combined.iter().map(|hit| hit.ctx_idx).collect();
    assert_eq!(contexts, vec![1, 0, 2, 0]);
    assert_eq!(combined[3].raw_score, 5);
}

/// The fields of NCBI's `BlastHSP` that `Blast_HSPListsMerge` reads or writes, for the
/// transcription of NCBI's merge below; `id` names the record.
#[derive(Clone, Debug, PartialEq)]
struct NcbiHsp {
    id: u32,
    context: i32,
    q_offset: i32,
    q_end: i32,
    s_offset: i32,
    s_end: i32,
    q_gapped_start: i32,
    s_gapped_start: i32,
    s_frame: i8,
    score: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353 (`ScoreCompareHSPs`)
fn ncbi_score_compare_hsps(hsp1: &NcbiHsp, hsp2: &NcbiHsp) -> std::cmp::Ordering {
    hsp2.score
        .cmp(&hsp1.score)
        .then(hsp1.s_offset.cmp(&hsp2.s_offset))
        .then(hsp2.s_end.cmp(&hsp1.s_end))
        .then(hsp1.q_offset.cmp(&hsp2.q_offset))
        .then(hsp2.q_end.cmp(&hsp1.q_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1488-1531 (`s_BlastMergeTwoHSPs`)
// A statement-by-statement transcription with `allow_gap` FALSE; `None` where `(int)` gets a
// value that is not finite (undefined in C).
fn ncbi_merge_two_hsps(hsp1: &mut NcbiHsp, hsp2: &NcbiHsp) -> Option<bool> {
    if hsp1.s_offset - hsp2.s_offset - hsp1.q_offset + hsp2.q_offset != 0 {
        return Some(false);
    }
    if hsp1.s_frame != hsp2.s_frame {
        return Some(false);
    }
    let contained =
        |a: i32, b: i32, c: i32, d: i32, e: i32, f: i32| (a <= c && b >= c) && (d <= f && e >= f);
    if contained(
        hsp1.q_offset,
        hsp1.q_end,
        hsp2.q_offset,
        hsp1.s_offset,
        hsp1.s_end,
        hsp2.s_offset,
    ) || contained(
        hsp1.q_offset,
        hsp1.q_end,
        hsp2.q_end,
        hsp1.s_offset,
        hsp1.s_end,
        hsp2.s_end,
    ) {
        let score_density = f64::from(hsp1.score + hsp2.score) * 1.0
            / f64::from((hsp1.q_end - hsp1.q_offset) + (hsp2.q_end - hsp2.q_offset));
        hsp1.q_offset = hsp1.q_offset.min(hsp2.q_offset);
        hsp1.s_offset = hsp1.s_offset.min(hsp2.s_offset);
        hsp1.q_end = hsp1.q_end.max(hsp2.q_end);
        hsp1.s_end = hsp1.s_end.max(hsp2.s_end);
        if hsp2.score > hsp1.score {
            hsp1.q_gapped_start = hsp2.q_gapped_start;
            hsp1.s_gapped_start = hsp2.s_gapped_start;
            hsp1.score = hsp2.score;
        }
        let product = score_density * f64::from(hsp1.q_end - hsp1.q_offset);
        if !product.is_finite() {
            return None;
        }
        hsp1.score = (product as i32).max(hsp1.score);
        return Some(true);
    }
    Some(false)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2876-2911,2961-3002,3029
// and :2758-2766, :2225-2247 (`Blast_HSPListPurgeNullHSPs`), :1355-1383
// `Blast_HSPListsMerge` of a split subject (`contexts_per_query` < 0), its pointer arrays as
// indices into `arena`, statement by statement; the qsort is glibc's merge sort (stable), the
// sort of `slice::sort_by`. `partition` false leaves out the swaps (the order LOSAT merged in
// before). `None`: a merge NCBI leaves undefined.
fn ncbi_hsp_lists_merge(
    arena: &mut [NcbiHsp],
    combined: Vec<usize>,
    new: Vec<usize>,
    split_offset: i32,
    chunk_overlap_size: i32,
    partition: bool,
) -> Option<Vec<usize>> {
    let mut combined_hsp_array = combined;
    let mut hsp_array: Vec<Option<usize>> = new.into_iter().map(Some).collect();
    let (mut hspcnt1, mut hspcnt2) = (0usize, 0usize);
    for index1 in 0..combined_hsp_array.len() {
        let hsp1 = combined_hsp_array[index1];
        if arena[hsp1].s_end > split_offset {
            if partition {
                let hsp_var = combined_hsp_array[hspcnt1];
                combined_hsp_array[hspcnt1] = hsp1;
                combined_hsp_array[index1] = hsp_var;
            }
            hspcnt1 += 1;
        }
    }
    for index2 in 0..hsp_array.len() {
        let hsp2 = hsp_array[index2].expect("no deleted HSP yet");
        if arena[hsp2].s_offset < split_offset + chunk_overlap_size {
            if partition {
                let hsp_var = hsp_array[hspcnt2];
                hsp_array[hspcnt2] = Some(hsp2);
                hsp_array[index2] = hsp_var;
            }
            hspcnt2 += 1;
        }
    }
    // Without the swaps, the strip HSPs are those that pass the test, in list order.
    let strip1: Vec<usize> = (0..combined_hsp_array.len())
        .filter(|&i| partition || arena[combined_hsp_array[i]].s_end > split_offset)
        .take(hspcnt1)
        .collect();
    let strip2: Vec<usize> = (0..hsp_array.len())
        .filter(|&i| {
            partition
                || arena[hsp_array[i].expect("not deleted")].s_offset
                    < split_offset + chunk_overlap_size
        })
        .take(hspcnt2)
        .collect();
    if hspcnt1 > 0 && hspcnt2 > 0 {
        for &index1 in &strip1 {
            let hsp1 = combined_hsp_array[index1];
            for &index2 in &strip2 {
                let Some(hsp2) = hsp_array[index2] else {
                    continue;
                };
                if arena[hsp1].context != arena[hsp2].context {
                    continue;
                }
                let end_diag = arena[hsp1].q_end - arena[hsp1].s_end;
                let start_diag = arena[hsp2].q_offset - arena[hsp2].s_offset;
                if (end_diag - start_diag).abs() < 10 {
                    let second = arena[hsp2].clone();
                    if ncbi_merge_two_hsps(&mut arena[hsp1], &second)? {
                        hsp_array[index2] = None;
                    }
                }
            }
        }
        let mut index1 = 0usize;
        for index in 0..hsp_array.len() {
            if hsp_array[index].is_some() {
                hsp_array[index1] = hsp_array[index];
                index1 += 1;
            }
        }
        hsp_array.truncate(index1);
    }
    for hsp in hsp_array.into_iter().flatten() {
        combined_hsp_array.push(hsp);
    }
    let sorted = combined_hsp_array.windows(2).all(|pair| {
        ncbi_score_compare_hsps(&arena[pair[0]], &arena[pair[1]]) != std::cmp::Ordering::Greater
    });
    if !sorted {
        combined_hsp_array.sort_by(|&a, &b| ncbi_score_compare_hsps(&arena[a], &arena[b]));
    }
    Some(combined_hsp_array)
}

fn ungapped_from_ncbi(hsp: &NcbiHsp) -> UngappedHit {
    let mut hit = make_int4_chunk_hit(
        hsp.context as usize,
        hsp.q_offset as usize,
        hsp.q_end,
        hsp.s_offset as usize,
        hsp.s_end,
        hsp.score,
    );
    hit.q_seed_off = hsp.q_gapped_start as u32;
    hit.s_seed_off = hsp.s_gapped_start as u32;
    hit.s_frame = hsp.s_frame;
    hit.hsp_list_order = hsp.id;
    hit
}

fn ncbi_from_ungapped(hit: &UngappedHit) -> NcbiHsp {
    NcbiHsp {
        id: hit.hsp_list_order,
        context: hit.ctx_idx as i32,
        q_offset: hit.q_aa_start as i32,
        q_end: hit.q_end_int4(),
        s_offset: hit.s_aa_start as i32,
        s_end: hit.s_end_int4(),
        q_gapped_start: hit.q_seed_off as i32,
        s_gapped_start: hit.s_seed_off as i32,
        s_frame: hit.s_frame,
        score: hit.raw_score,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-3031 (`Blast_HSPListsMerge`)
// The port against the transcription of NCBI's merge above, on random lists of two chunks
// (split 200, overlap 100) sorted by score as NCBI's are: HSPs of length -2 to 12, groups of
// records that tie on every key of the score sort (the same coordinates and score in up to
// four contexts), HSPs in and out of the overlap strip on both sides, HSPs of the later chunk
// on the diagonal of one of the earlier chunk (merged and deleted). Every list the
// transcription gives is the port's, record for record in its order; where the transcription
// meets an undefined merge, the port refuses. Some lists come out in another order than
// without the swaps (the order LOSAT merged in before), some merges delete HSPs and some are
// refused.
#[test]
fn chunk_merge_is_ncbis_merge() {
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move |n: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % n
    };
    let (split, overlap) = (200i32, 100i32);
    let (mut reordered, mut merged, mut refused, mut compared) = (0usize, 0usize, 0usize, 0usize);
    for case in 0..3000 {
        let mut arena: Vec<NcbiHsp> = Vec::new();
        let push_group = |arena: &mut Vec<NcbiHsp>,
                          next: &mut dyn FnMut(u64) -> u64,
                          q_offset: i32,
                          s_offset: i32,
                          length: i32,
                          s_frame: i8|
         -> Vec<usize> {
            let score = 1 + next(12) as i32;
            let copies = 1 + next(4) as usize;
            let first_context = next(3) as i32;
            (0..copies)
                .map(|copy| {
                    arena.push(NcbiHsp {
                        id: arena.len() as u32,
                        context: first_context + copy as i32,
                        q_offset,
                        q_end: q_offset + length,
                        s_offset,
                        s_end: s_offset + length,
                        q_gapped_start: q_offset + length.max(1) / 2,
                        s_gapped_start: s_offset + length.max(1) / 2,
                        s_frame,
                        score,
                    });
                    arena.len() - 1
                })
                .collect()
        };
        let length = |next: &mut dyn FnMut(u64) -> u64| next(15) as i32 - 2;
        let frame = |next: &mut dyn FnMut(u64) -> u64| if next(4) == 0 { -1 } else { 1 };
        let mut first: Vec<usize> = Vec::new();
        for _ in 0..1 + next(6) {
            // A quarter of the groups start at subject offset 0 or 1, so that their subject
            // ends (offset + length, length -2 and up) can be below 0: an `Int4` end and a
            // `u32` end tell apart in the strip test and in the final sort.
            let s_offset = if next(4) == 0 {
                next(2) as i32
            } else {
                170 + next(130) as i32
            };
            let q_offset = (s_offset - 170 + next(4) as i32).max(0);
            let (len, s_frame) = (length(&mut next), frame(&mut next));
            first.extend(push_group(
                &mut arena, &mut next, q_offset, s_offset, len, s_frame,
            ));
        }
        let mut second: Vec<usize> = Vec::new();
        for _ in 0..1 + next(6) {
            let (q_offset, s_offset, s_frame) = if next(2) == 0 {
                // On the diagonal of an HSP of the earlier chunk, from its start on.
                let base = arena[first[next(first.len() as u64) as usize]].clone();
                let shift = next(4) as i32;
                (base.q_offset + shift, base.s_offset + shift, base.s_frame)
            } else {
                let s_offset = 200 + next(130) as i32;
                (s_offset - 170 + next(4) as i32, s_offset, frame(&mut next))
            };
            let len = length(&mut next);
            second.extend(push_group(
                &mut arena, &mut next, q_offset, s_offset, len, s_frame,
            ));
        }
        let by_score = |arena: &[NcbiHsp], list: &mut Vec<usize>| {
            list.sort_by(|&a, &b| ncbi_score_compare_hsps(&arena[a], &arena[b]));
        };
        by_score(&arena, &mut first);
        by_score(&arena, &mut second);

        let mut ncbi_arena = arena.clone();
        let expected = ncbi_hsp_lists_merge(
            &mut ncbi_arena,
            first.clone(),
            second.clone(),
            split,
            overlap,
            true,
        );
        let mut combined: Vec<UngappedHit> = first
            .iter()
            .map(|&i| ungapped_from_ncbi(&arena[i]))
            .collect();
        let incoming: Vec<UngappedHit> = second
            .iter()
            .map(|&i| ungapped_from_ncbi(&arena[i]))
            .collect();
        let result = merge_tblastx_subject_chunk_hits(&mut combined, incoming, split as usize, 100);
        let Some(expected) = expected else {
            assert!(result.is_err(), "case {case}: NCBI's merge is undefined");
            refused += 1;
            continue;
        };
        assert!(result.is_ok(), "case {case}: NCBI's merge is defined");
        let actual: Vec<NcbiHsp> = combined.iter().map(ncbi_from_ungapped).collect();
        let expected: Vec<NcbiHsp> = expected.iter().map(|&i| ncbi_arena[i].clone()).collect();
        assert_eq!(actual, expected, "case {case}");
        compared += 1;
        merged += usize::from(expected.len() < first.len() + second.len());
        let mut unswapped_arena = arena.clone();
        let unswapped =
            ncbi_hsp_lists_merge(&mut unswapped_arena, first, second, split, overlap, false)
                .expect("the same merges");
        let unswapped: Vec<u32> = unswapped.iter().map(|&i| unswapped_arena[i].id).collect();
        reordered += usize::from(unswapped != expected.iter().map(|h| h.id).collect::<Vec<_>>());
    }
    assert!(compared > 1000 && merged > 100 && refused > 10 && reordered > 10);
}

/// The TBLASTX engine of a strategy that another application runs (R2-B1, review 2 of F2).
fn run_as(running: &str, options: &[&str]) -> TblastxArgs {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/strategy/import/ss/tblastx.db.ss");
    let mut argv = vec![
        "losat",
        running,
        "-import_search_strategy",
        path.to_str().unwrap(),
    ];
    argv.extend_from_slice(options);
    let cli = crate::cli::try_parse_from(argv).unwrap();
    let command =
        crate::report::archive::read::import::apply(cli.command, None, &mut Vec::new()).unwrap();
    let crate::cli::Commands::Tblastx(args) = command else {
        panic!("{running} running a tblastx strategy");
    };
    args
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2947-2989
// The formatting handler's warnings of blastn (and of blastx) running a tblastx strategy carry the
// running application's name.
#[test]
fn the_warnings_of_another_applications_run_carry_its_name() {
    for running in ["blastn", "blastx"] {
        let args = run_as(running, &["-max_target_seqs", "2"]);
        assert_eq!(args.app().ncbi_name(), running);
        let mut posted = Vec::new();
        check_ncbi_options(&args, &cli_outfmts(&args).unwrap(), &mut posted).unwrap();
        let posted = String::from_utf8(posted).unwrap();
        assert!(
            posted.contains(&format!(
                "Warning: [{running}] Examining 5 or more matches is recommended\n"
            )),
            "{running}: {posted}"
        );
        assert!(!posted.contains("[tblastx]"), "{running}: {posted}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2875-2878
// The SAM check is the running application's: blastx stops with NCBI's error, blastn passes it and
// writes the SAM report of the TBLASTX strategy (`tblastx::sam`).
#[test]
fn the_sam_check_of_another_applications_run_is_the_running_applications() {
    let checked = |running: &str| {
        let args = run_as(running, &["-outfmt", "17"]);
        check_ncbi_options(&args, &cli_outfmts(&args).unwrap(), &mut Vec::new())
            .map(|_| ())
            .map_err(|error| error.to_string())
    };
    assert!(checked("blastx")
        .unwrap_err()
        .contains("SAM format is only applicable to blastn"));
    assert_eq!(checked("blastn"), Ok(()));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:42-65,137 (the table of a window
// below 0 is the one-hit one, its offset starts at the window and its cells at 0) and
// c++/src/algo/blast/core/aa_ungapped.c:775-783 (a word is extended when `diff >= 0`)
#[test]
fn a_fresh_table_of_a_negative_window_skips_the_words_before_the_window() {
    for (window, first) in [(-1, 1), (-40, 40)] {
        let offset = word_finder::WideOffset::new(window);
        let cell = DiagStruct::default();
        assert_eq!(offset.get(), i64::from(window));
        for s_off in [0, first - 1] {
            assert!(
                !word_finder::WideOffset::admits(offset.get(), cell.last_hit(), s_off),
                "window {window}, word at {s_off}"
            );
        }
        for s_off in [first, first + 1, 5_000] {
            assert!(
                word_finder::WideOffset::admits(offset.get(), cell.last_hit(), s_off),
                "window {window}, word at {s_off}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-185 (`Blast_ExtendWordExit`
// adds the subject's length and the window) and c++/src/algo/blast/core/aa_ungapped.c:775-800
#[test]
fn a_table_of_a_negative_window_carries_its_cells_to_the_next_call() {
    let window = -40;
    let wordsize = 3;
    let mut cells = vec![DiagStruct::default(); 4];
    let offset = word_finder::WideOffset::new(window);
    // An extension that ended at 300 stores `300 - (wordsize - 1) + offset`.
    let stored = word_finder::WideOffset::last_hit(offset.get(), 300, wordsize);
    assert_eq!(stored, 258);
    cells[0].set_last_hit(stored);
    // A long call moves the offset past the cell: nothing it holds blocks a word of the next.
    let after_long = word_finder::WideOffset::new(window);
    after_long.exit(window, &mut cells, 1_000);
    assert_eq!(after_long.get(), -40 + 1_000 - 40);
    assert!(word_finder::WideOffset::admits(after_long.get(), 258, 0));
    // A short call moves it by less than the cell's lead: the next call's words before 328 on
    // that diagonal are blocked, and those before 70 on a cell that was never written: a table
    // made fresh for the call would block neither.
    offset.exit(window, &mut cells, 10);
    assert_eq!(offset.get(), -70);
    assert!(!word_finder::WideOffset::admits(offset.get(), 258, 327));
    assert!(word_finder::WideOffset::admits(offset.get(), 258, 328));
    let unwritten = DiagStruct::default().last_hit();
    assert!(!word_finder::WideOffset::admits(
        offset.get(),
        unwritten,
        69
    ));
    assert!(word_finder::WideOffset::admits(offset.get(), unwritten, 70));
    // The cells themselves stay as they were stored.
    assert_eq!(cells[0].last_hit(), 258);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:90-106,162-173 (the table is
// cleared to `-window` and the offset goes back to the window once the offset reaches
// a quarter of `INT4_MAX`)
#[test]
fn a_table_of_a_negative_window_is_cleared_when_its_offset_reaches_a_quarter_of_int4_max() {
    let window = -1;
    let mut cells = vec![DiagStruct::default(); 4];
    cells[2].set_last_hit(123);
    cells[2].set_flag(1);
    let offset = word_finder::WideOffset::new(window);
    let mut calls: i32 = 0;
    while offset.get() < i64::from(i32::MAX / 4) {
        offset.exit(window, &mut cells, 5_000_000);
        calls += 1;
        assert!(calls < 200, "the offset rises by 4,999,999 a call");
    }
    assert_eq!(offset.get(), -1 + i64::from(calls) * 4_999_999);
    assert_eq!(
        cells[2].last_hit(),
        123,
        "no clear before the exit that finds the offset"
    );
    offset.exit(window, &mut cells, 5_000_000);
    assert_eq!(offset.get(), -1);
    for cell in &cells {
        assert_eq!((cell.last_hit(), cell.flag()), (1, 0));
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-61 (a 31-bit signed
// `last_hit`, truncated when stored)
#[test]
fn a_table_of_a_negative_window_stores_the_31_low_bits_of_the_cell() {
    let mut cell = DiagStruct::default();
    cell.set_last_hit(word_finder::WideOffset::last_hit(1 << 30, 10, 3));
    assert_eq!(cell.last_hit(), -1_073_741_816);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-185 and
// c++/src/algo/blast/core/aa_ungapped.c:775-800: sums that overflow `Int4` there (undefined) are
// exact here, and nothing panics for any window
#[test]
fn a_table_of_the_smallest_window_runs_without_overflow() {
    let window = i32::MIN;
    let mut cells = vec![DiagStruct::default(); 1];
    let offset = word_finder::WideOffset::new(window);
    assert_eq!(offset.get(), i64::from(i32::MIN));
    // `last_hit - diag_offset` is 2^31: no word of a subject of fewer than 2^31 letters is extended.
    for s_off in [0, 1, 5_000_000] {
        assert!(!word_finder::WideOffset::admits(offset.get(), 0, s_off));
    }
    offset.exit(window, &mut cells, 1_000);
    assert_eq!(offset.get(), 2 * i64::from(i32::MIN) + 1_000);
    offset.exit(window, &mut cells, 1_000);
    assert_eq!(offset.get(), 3 * i64::from(i32::MIN) + 2_000);
    assert!(!word_finder::WideOffset::admits(offset.get(), 0, 5_000_000));
    // The cell keeps the 31 low bits of the exact sum.
    let mut cell = DiagStruct::default();
    cell.set_last_hit(word_finder::WideOffset::last_hit(
        -(1 << 32) + 1_000,
        100,
        3,
    ));
    assert_eq!(cell.last_hit(), 1_098);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:775-783 (`diff` is a `Uint4`
// minus an `Int4` read as an `Int4`): where NCBI's arithmetic is defined, the exact offset gives
// the test of the one-hit scan of the window 0 (`scan_chunk_one_hit`)
#[test]
fn the_offset_test_of_a_negative_window_is_the_one_hit_test_where_int4_holds() {
    let narrow = |offset: i32, last_hit: i32, s_off: u32| {
        s_off.wrapping_sub(last_hit.wrapping_sub(offset) as u32) as i32 >= 0
    };
    let offsets: [i32; 10] = [
        -2_000_000_000,
        -1_073_741_000,
        -1_000,
        -41,
        -1,
        0,
        7,
        1_000,
        536_870_911,
        1_000_000_000,
    ];
    let last_hits: [i32; 6] = [0, 5, -5, 258, (1 << 30) - 1, -(1 << 30)];
    let words = [0u32, 1, 39, 40, 100, 328, 1 << 20, (1 << 30) + 9];
    for offset in offsets {
        for last_hit in last_hits {
            if last_hit.checked_sub(offset).is_none() {
                continue;
            }
            for s_off in words {
                assert_eq!(
                    word_finder::WideOffset::admits(i64::from(offset), last_hit, s_off),
                    narrow(offset, last_hit, s_off),
                    "offset {offset}, cell {last_hit}, word {s_off}"
                );
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
// A negative chunk length moves `diag_offset` back by the chunk length plus window; at or above
// INT4_MAX / 4 the table clears to `DiagStruct::clear(window)` and `diag_offset` resets to window.
#[test]
fn advance_tblastx_diag_offset_handles_negative_length_and_clear() {
    let window = 40;
    let mut diag_array = vec![DiagStruct::default(); 4];
    diag_array[0].set_last_hit(123);
    diag_array[0].set_flag(1);
    let mut diag_offset = 40;

    advance_tblastx_diag_offset(&mut diag_offset, &mut diag_array, window, -100);
    assert_eq!(diag_offset, -20);
    assert_eq!((diag_array[0].last_hit(), diag_array[0].flag()), (123, 1));

    let clear_at = i32::MAX / 4;
    for (start, len, want_offset, cleared) in [
        (clear_at - 1, -5, clear_at - 1 - 5 + window, false),
        (clear_at, -5, window, true),
    ] {
        diag_offset = start;
        advance_tblastx_diag_offset(&mut diag_offset, &mut diag_array, window, len);
        assert_eq!(diag_offset, want_offset);
        if cleared {
            assert!(diag_array.iter().all(|d| *d == DiagStruct::clear(window)));
        } else {
            assert_eq!(diag_array[0].last_hit(), 123);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-185
// A negative window over overlapping hard masks (negative signed chunk length) moves `WideOffset`
// back by the chunk length plus window without resetting or losing the carried table, and a
// >15 Mbp subject (`s_len / 3 > MAX_DBSEQ_LEN`) keeps that carried table across its chunks.
#[test]
fn wide_offset_exit_handles_negative_chunk_length_and_long_subject_carry() {
    let window = -40;
    let mut cells = vec![DiagStruct::default(); 4];
    let offset = word_finder::WideOffset::new(window);
    assert_eq!(offset.get(), -40);
    // Overlapping hard-mask chunk of signed length -25 moves offset by -25 + (-40).
    offset.exit(window, &mut cells, -25);
    assert_eq!(offset.get(), -105);
    // A 5,000,000-aa chunk of a >15 Mbp subject advances the carried table as in NCBI.
    offset.exit(window, &mut cells, 5_000_000);
    assert_eq!(offset.get(), -105 + 5_000_000 - 40);
}

// Line RG series D, D-SO-5: the held writes after the loop (`BatchReports::write_batches`) take
// the batches as the stream takes them: each batch's hits, query statistics and flags, its hits
// numbered from its first query (the loop numbers the held hits from the run's first query,
// `run_batch_loop`), for batches with and without hits, of one query and of several.
#[test]
fn held_results_are_split_into_the_batches_the_stream_takes() {
    fn hsp(q_idx: u32, s_idx: u32, raw_score: i32) -> TblastxHsp {
        TblastxHsp {
            hit: crate::common::Hit {
                identity: 100.0,
                length: 10,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 30,
                s_start: 1,
                s_end: 30,
                e_value: 1e-9,
                bit_score: 50.0,
                num_ident: 10,
                query_frame: 1,
                query_length: 0,
                q_idx,
                s_idx,
                raw_score,
                sort_query_offset: 0,
                sort_query_end: 10,
                sort_subject_offset: 0,
                sort_subject_end: 10,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 10,
            },
            subject_frame: 1,
            num: 1,
        }
    }
    let key = |hsp: &TblastxHsp| (hsp.hit.q_idx, hsp.hit.s_idx, hsp.hit.raw_score);
    // The batches as the search returns them: queries counted from 0 in each batch.
    let mut state = 0x2545_f491_u32;
    let mut next = |bound: u32| {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        state % bound
    };
    for counts in [vec![1, 1], vec![3, 1, 2], vec![2, 4, 1, 3], vec![5, 5]] {
        let batches: Vec<(Vec<TblastxHsp>, Vec<bool>)> = counts
            .iter()
            .map(|&count| {
                let mut hits = Vec::new();
                for query in 0..count {
                    for _ in 0..next(4) {
                        hits.push(hsp(query, next(50), next(1000) as i32));
                    }
                }
                let searched = (0..count).map(|_| next(2) == 0).collect();
                (hits, searched)
            })
            .collect();
        // What the loop holds: every batch's hits numbered from the run's first query.
        let mut done = SearchedBatches::default();
        let mut start = 0;
        for (hits, searched) in &batches {
            let count = searched.len();
            done.hits.extend(hits.iter().cloned().map(|mut hsp| {
                hsp.hit.q_idx += start as u32;
                hsp
            }));
            done.query_stats
                .extend((0..count).map(|query| TblastxQueryStats {
                    eff_searchsp: (start + query) as i64,
                    ..TblastxQueryStats::default()
                }));
            done.searched.extend(searched);
            done.spans.push(start..start + count);
            start += count;
        }
        let split = done.into_batches();
        assert_eq!(split.len(), batches.len(), "{counts:?}");
        let mut start = 0;
        for ((hits, searched), batch) in batches.iter().zip(&split) {
            let count = searched.len();
            assert_eq!(
                batch.hits.iter().map(key).collect::<Vec<_>>(),
                hits.iter().map(key).collect::<Vec<_>>(),
                "{counts:?}"
            );
            assert_eq!(&batch.searched, searched, "{counts:?}");
            assert_eq!(
                batch
                    .query_stats
                    .iter()
                    .map(|stats| stats.eff_searchsp)
                    .collect::<Vec<_>>(),
                (start..start + count)
                    .map(|query| query as i64)
                    .collect::<Vec<_>>(),
                "{counts:?}"
            );
            assert_eq!(batch.spans, vec![0..count], "{counts:?}");
            start += count;
        }
    }
    assert!(SearchedBatches::default().into_batches().is_empty());
}
