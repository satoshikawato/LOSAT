use super::*;
use crate::stats::KarlinParams;

fn make_context() -> QueryContext {
    QueryContext {
        q_idx: 0,
        frame: 1,
        aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
        aa_seq_nomask: None,
        aa_len: 5,
        orig_len: 15,
        frame_base: 0,
        is_valid: true,
        karlin_params: KarlinParams::default(),
    }
}

fn make_subject_frame() -> QueryFrame {
    QueryFrame {
        frame: 1,
        aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
        aa_seq_nomask: None,
        aa_len: 5,
        orig_len: 15,
        seg_masks: Vec::new(),
    }
}

fn make_init_hsp(q_start: i32, q_end: i32, s_start: i32, s_end: i32, score: i32) -> InitHSP {
    InitHSP {
        q_start_absolute: q_start,
        q_end_absolute: q_end,
        s_start,
        s_end,
        q_seed_absolute: q_start,
        s_seed: s_start,
        score,
        ctx_idx: 0,
        s_f_idx: 0,
        q_idx: 0,
        s_idx: 0,
        q_frame: 1,
        s_frame: 1,
        q_orig_len: 15,
        s_orig_len: 15,
    }
}

fn make_ungapped_for_sort(raw_score: i32, hsp_list_order: usize) -> UngappedHit {
    UngappedHit {
        q_idx: hsp_list_order as u32,
        s_idx: 0,
        ctx_idx: 0,
        s_f_idx: 0,
        q_frame: 1,
        s_frame: 1,
        q_aa_start: 10,
        q_aa_end: 20,
        s_aa_start: 30,
        s_aa_end: 40,
        q_seed_off: 10,
        s_seed_off: 30,
        q_orig_len: 100,
        s_orig_len: 100,
        raw_score,
        e_value: 0.0,
        num_ident: 0,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
        hsp_list_order: u32::try_from(hsp_list_order).expect("hsp_list_order fits u32"),
        ordering_method: 0,
        linked_set: false,
        start_of_chain: false,
        link_id: 0,
        chain_next_link_id: None,
        hsp_link_num: 0,
        num: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-310
// score_compare_match breaks score/s_start/length ties on the ABSOLUTE
// q_start, so equal context-relative HSPs of different query frames leave
// the word finder in context order whatever their scan order was.
#[test]
fn test_sort_init_hsps_orders_equal_relative_hits_by_absolute_query_start() {
    let mut ctx2_first = make_init_hsp(20 + 7, 20 + 17, 519, 529, 42);
    ctx2_first.ctx_idx = 2;
    let mut ctx0_second = make_init_hsp(7, 17, 519, 529, 42);
    ctx0_second.ctx_idx = 0;
    let mut hsps = vec![ctx2_first, ctx0_second];
    sort_init_hsps_by_score_ncbi(&mut hsps);
    assert_eq!(
        hsps.iter().map(|h| h.ctx_idx).collect::<Vec<_>>(),
        vec![0, 2]
    );
}

#[test]
fn test_get_ungapped_hsp_list_sorts_like_ncbi_score_compare_match() {
    let contexts = vec![make_context()];
    let subject_frames = vec![make_subject_frame()];
    let hits = get_ungapped_hsp_list(
        vec![
            make_init_hsp(4, 10, 30, 36, 100),
            make_init_hsp(2, 10, 20, 28, 100),
            make_init_hsp(1, 8, 20, 27, 100),
            make_init_hsp(3, 11, 20, 28, 120),
        ],
        &contexts,
        &subject_frames,
    );

    let ordered: Vec<(i32, usize, usize, usize)> = hits
        .iter()
        .map(|hit| {
            (
                hit.raw_score,
                hit.s_aa_start as usize,
                (hit.q_aa_end as usize).saturating_sub(hit.q_aa_start as usize),
                hit.q_aa_start as usize,
            )
        })
        .collect();

    assert_eq!(
        ordered,
        vec![
            (120, 20, 8, 3),
            (100, 20, 8, 2),
            (100, 20, 7, 1),
            (100, 30, 6, 4)
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379
#[test]
fn test_get_ungapped_hsp_list_preserves_insertion_order_when_score_sorted() {
    let contexts = vec![
        QueryContext {
            q_idx: 0,
            frame: 1,
            aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
            aa_seq_nomask: None,
            aa_len: 5,
            orig_len: 15,
            frame_base: 0,
            is_valid: true,
            karlin_params: KarlinParams::default(),
        },
        QueryContext {
            q_idx: 1,
            frame: 1,
            aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
            aa_seq_nomask: None,
            aa_len: 5,
            orig_len: 15,
            frame_base: 100,
            is_valid: true,
            karlin_params: KarlinParams::default(),
        },
    ];
    let subject_frames = vec![make_subject_frame()];
    let hits = get_ungapped_hsp_list(
        vec![
            InitHSP {
                q_start_absolute: 101,
                q_end_absolute: 106,
                s_start: 20,
                s_end: 25,
                q_seed_absolute: 101,
                s_seed: 20,
                score: 100,
                ctx_idx: 1,
                s_f_idx: 0,
                q_idx: 1,
                s_idx: 0,
                q_frame: 1,
                s_frame: 1,
                q_orig_len: 15,
                s_orig_len: 15,
            },
            InitHSP {
                q_start_absolute: 1,
                q_end_absolute: 6,
                s_start: 20,
                s_end: 25,
                q_seed_absolute: 1,
                s_seed: 20,
                score: 100,
                ctx_idx: 0,
                s_f_idx: 0,
                q_idx: 0,
                s_idx: 0,
                q_frame: 1,
                s_frame: 1,
                q_orig_len: 15,
                s_orig_len: 15,
            },
        ],
        &contexts,
        &subject_frames,
    );

    let ordered_qidx: Vec<u32> = hits.iter().map(|hit| hit.q_idx).collect();
    assert_eq!(ordered_qidx, vec![1, 0]);
}

#[test]
fn test_ncbi_score_sort_preserves_current_hsp_list_order_for_equal_hits() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1355
    let mut hits = vec![
        make_ungapped_for_sort(90, 0),
        make_ungapped_for_sort(100, 2),
        make_ungapped_for_sort(100, 1),
    ];

    sort_ungapped_hits_by_score_ncbi(&mut hits);

    let ordered: Vec<(i32, usize)> = hits
        .iter()
        .map(|hit| (hit.raw_score, hit.hsp_list_order as usize))
        .collect();
    assert_eq!(ordered, vec![(100, 2), (100, 1), (90, 0)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2392
#[test]
#[should_panic(expected = "NCBI BLAST requires context-relative ungapped q_start >= 0")]
fn test_get_ungapped_hsp_list_fails_fast_on_negative_context_relative_query_start() {
    let contexts = vec![QueryContext {
        q_idx: 0,
        frame: 1,
        aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
        aa_seq_nomask: None,
        aa_len: 5,
        orig_len: 15,
        frame_base: 10,
        is_valid: true,
        karlin_params: KarlinParams::default(),
    }];
    let subject_frames = vec![make_subject_frame()];
    let invalid = InitHSP {
        q_start_absolute: 5,
        q_end_absolute: 10,
        s_start: 20,
        s_end: 25,
        q_seed_absolute: 5,
        s_seed: 20,
        score: 100,
        ctx_idx: 0,
        s_f_idx: 0,
        q_idx: 0,
        s_idx: 0,
        q_frame: 1,
        s_frame: 1,
        q_orig_len: 15,
        s_orig_len: 15,
    };

    let _ = get_ungapped_hsp_list(vec![invalid], &contexts, &subject_frames);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4730-4732
// The list passed in by the caller (LOSAT_X_SEEDCOMPACT reserves one hit per init HSP) gives
// the list of `get_ungapped_hsp_list`, every field (Debug, so E-values by their shortest
// round-trip text), in the caller's allocation; an empty input keeps the reservation.
#[test]
fn test_get_ungapped_hsp_list_into_a_reserved_list_equals_get_ungapped_hsp_list() {
    let contexts = vec![make_context()];
    let subject_frames = vec![make_subject_frame()];
    let init_hsps = vec![
        make_init_hsp(4, 10, 30, 36, 100),
        make_init_hsp(2, 10, 20, 28, 100),
        make_init_hsp(1, 8, 20, 27, 100),
        make_init_hsp(3, 11, 20, 28, 120),
    ];
    let expected = get_ungapped_hsp_list(init_hsps.clone(), &contexts, &subject_frames);
    let reserved = Vec::with_capacity(init_hsps.len());
    let buffer = reserved.as_ptr();
    let hits = get_ungapped_hsp_list_into(init_hsps, &contexts, &subject_frames, reserved);
    assert_eq!(format!("{hits:?}"), format!("{expected:?}"));
    assert_eq!((hits.len(), hits.capacity()), (4, 4));
    assert_eq!(hits.as_ptr(), buffer, "the caller's allocation");
    let empty = get_ungapped_hsp_list_into(
        Vec::new(),
        &contexts,
        &subject_frames,
        Vec::with_capacity(3),
    );
    assert!(empty.is_empty() && empty.capacity() == 3);
}

// No NCBI counterpart: the caller-given list must be empty (`hsp_list_order` counts from its
// length).
#[test]
#[should_panic(expected = "get_ungapped_hsp_list_into: the list passed in must be empty")]
fn test_get_ungapped_hsp_list_into_rejects_a_list_with_hits() {
    let contexts = vec![make_context()];
    let subject_frames = vec![make_subject_frame()];
    let one = get_ungapped_hsp_list(
        vec![make_init_hsp(1, 8, 20, 27, 100)],
        &contexts,
        &subject_frames,
    );
    let _ = get_ungapped_hsp_list_into(
        vec![make_init_hsp(2, 10, 20, 28, 100)],
        &contexts,
        &subject_frames,
        one,
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4760-4772
// The one-hit conversion keeps the ends `start + length` of HSPs of length 0, -1 and -2 (one
// ending at -2, below the frame's start) and sorts with the `Int4` comparator: of two HSPs of
// score 9 that start at subject 0, the one ending at 3 comes before the one ending at -2
// (`subject.end` descending). The two-hit lists' comparator reads the subject's offset and end as
// `Int4` too (a chunk's HSP before the chunk, M4b-T2), so it orders them the same way. The HSPs of
// positive length are `get_ungapped_hsp_list`'s, field for field.
#[test]
fn int4_conversion_keeps_ends_at_or_before_the_start_and_sorts_them_as_int4() {
    let mut context = make_context();
    context.frame_base = 100;
    let contexts = vec![context];
    let subject_frames = vec![make_subject_frame()];
    let init_hsps = vec![
        make_init_hsp(100, 98, 0, -2, 9),
        make_init_hsp(101, 104, 0, 3, 9),
        make_init_hsp(103, 103, 2, 2, 4),
        make_init_hsp(102, 101, 1, 0, 12),
    ];
    let hits = get_ungapped_hsp_list_int4(init_hsps.clone(), &contexts, &subject_frames);
    let ends: Vec<(u32, i32, u32, i32, i32)> = hits
        .iter()
        .map(|hit| {
            (
                hit.q_aa_start,
                hit.q_end_int4(),
                hit.s_aa_start,
                hit.s_end_int4(),
                hit.raw_score,
            )
        })
        .collect();
    assert_eq!(
        ends,
        vec![
            (2, 1, 1, 0, 12),
            (1, 4, 0, 3, 9),
            (0, -2, 0, -2, 9),
            (3, 3, 2, 2, 4),
        ]
    );
    assert_eq!(
        hits.iter()
            .map(UngappedHit::is_int4_inverted)
            .collect::<Vec<_>>(),
        vec![true, false, true, true]
    );
    let mut by_u32 = hits.clone();
    sort_ungapped_hits_by_index_replay(&mut by_u32, score_compare_ungapped_hits_ncbi);
    assert_eq!(by_u32[1].s_end_int4(), 3);
    assert_eq!(by_u32[2].s_end_int4(), -2);

    let positive = vec![
        make_init_hsp(101, 104, 0, 3, 9),
        make_init_hsp(105, 109, 4, 8, 9),
        make_init_hsp(102, 110, 1, 9, 20),
    ];
    let expected = get_ungapped_hsp_list(positive.clone(), &contexts, &subject_frames);
    let actual = get_ungapped_hsp_list_int4(positive, &contexts, &subject_frames);
    crate::algorithm::tblastx::sum_stats_linking::x_hsp_inplace::x_assert_same_hits(
        "positive lengths",
        &actual,
        &expected,
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:150-189,1330-1353
// A chunk's HSPs before the chunk (masks stored out of order, M4b-T2): the one-hit HSP of a word at
// offset -1 that starts at 0 (its seed at -1) and one from -2, then a two-hit HSP from -1. The
// conversion keeps the subject's offsets as NCBI's `Int4` and each sort reads them so: of HSPs of
// one score, the one that starts lower comes first.
#[test]
fn hsps_before_their_chunk_convert_and_sort_as_int4() {
    let contexts = vec![make_context()];
    let subject_frames = vec![make_subject_frame()];
    let mut seed_before = make_init_hsp(0, 3, 0, 3, 9);
    seed_before.s_seed = -1;
    let one_hit = get_ungapped_hsp_list_int4(
        vec![seed_before, make_init_hsp(1, 4, -2, 1, 9)],
        &contexts,
        &subject_frames,
    );
    let offsets: Vec<(i32, i32, i32)> = one_hit
        .iter()
        .map(|hit| (hit.s_start_int4(), hit.s_end_int4(), hit.s_seed_int4()))
        .collect();
    assert_eq!(offsets, vec![(-2, 1, -2), (0, 3, -1)]);
    let two_hit = get_ungapped_hsp_list(
        vec![make_init_hsp(0, 3, 0, 3, 9), make_init_hsp(2, 3, -1, 0, 9)],
        &contexts,
        &subject_frames,
    );
    let starts: Vec<i32> = two_hit.iter().map(UngappedHit::s_start_int4).collect();
    assert_eq!(starts, vec![-1, 0]);
}
