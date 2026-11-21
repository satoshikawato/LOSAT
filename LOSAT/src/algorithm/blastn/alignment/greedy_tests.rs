use super::*;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:380-387,537
/// The test runs the port of `BLAST_GreedyAlign` (and `BLAST_AffineGreedyAlign`) on random
/// inputs and compares the alignments.
///
/// On random pairs of related sequences, with and without traceback: the alignment on a
/// scratch that earlier, unrelated alignments have used equals the one on a new scratch,
/// which is what running tracebacks ahead of time on other scratches relies on.
#[test]
fn x_greedy_alignment_on_a_used_scratch_equals_a_new_scratch() {
    fn run(
        query: &[u8],
        subject: &[u8],
        q_off: usize,
        s_off: usize,
        scoring: (i32, i32, i32, i32),
        x_drop: i32,
        scratch: &mut GreedyAlignScratch,
        traceback: bool,
    ) -> Option<(
        i32,
        i32,
        i32,
        i32,
        i32,
        i32,
        i32,
        Option<(Vec<GapAlignOpType>, Vec<i32>)>,
    )> {
        let core = greedy_gapped_alignment_internal(
            ContextQuery {
                buf: query,
                start: 0,
                len: query.len() as i32,
            },
            SeqPointer {
                buf: subject,
                start: 0,
                len: subject.len() as i32,
            },
            subject.len(),
            q_off as i32,
            s_off as i32,
            scoring.0,
            scoring.1,
            scoring.2,
            scoring.3,
            x_drop,
            scratch,
            traceback,
        );
        core.map(|c| {
            (
                c.q_start,
                c.q_end,
                c.s_start,
                c.s_end,
                c.q_seed_start,
                c.s_seed_start,
                c.score,
                c.edit_script
                    .map(|e| (e.op_type[..e.size].to_vec(), e.num[..e.size].to_vec())),
            )
        })
    }

    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1500);
    let mut rng = TestRng(0x2545_F491_4F6C_DD1D);
    let mut used = GreedyAlignScratch::new();
    let mut aligned = 0usize;
    let mut with_gaps = 0usize;
    for case in 0..cases {
        // a core, a mutated copy of it, and unrelated flanks on both
        let core_len = 20 + rng.below(if case % 10 == 0 { 6000 } else { 600 }) as usize;
        let substitution = [0u64, 2, 5, 20, 60, 150][rng.below(6) as usize]; // per mille
        let indel = [0u64, 0, 2, 10, 30][rng.below(5) as usize]; // per mille
        let ambiguity = if case % 7 == 0 { 8 } else { 0 }; // per mille
        let random_base = |rng: &mut TestRng| rng.below(4) as u8;
        let core: Vec<u8> = (0..core_len).map(|_| random_base(&mut rng)).collect();
        let mut query = Vec::new();
        let mut subject = Vec::new();
        for _ in 0..rng.below(40) {
            query.push(random_base(&mut rng));
        }
        for _ in 0..rng.below(40) {
            subject.push(random_base(&mut rng));
        }
        // positions of an unchanged core base in both sequences
        let mut anchors: Vec<(usize, usize)> = Vec::new();
        for &base in &core {
            let roll = rng.below(1000);
            if roll < indel {
                // drop from the subject
                query.push(base);
            } else if roll < 2 * indel {
                // extra base in the subject
                subject.push(random_base(&mut rng));
                anchors.push((query.len(), subject.len()));
                query.push(base);
                subject.push(base);
            } else if roll < 2 * indel + substitution {
                query.push(base);
                subject.push((base + 1 + rng.below(3) as u8) & 3);
            } else if roll < 2 * indel + substitution + ambiguity {
                // an ambiguity code in the query (never matches)
                query.push(4 + rng.below(12) as u8);
                subject.push(base);
            } else {
                anchors.push((query.len(), subject.len()));
                query.push(base);
                subject.push(base);
            }
        }
        for _ in 0..rng.below(40) {
            query.push(random_base(&mut rng));
        }
        for _ in 0..rng.below(40) {
            subject.push(random_base(&mut rng));
        }
        if anchors.is_empty() {
            continue;
        }
        let (q_off, s_off) = anchors[rng.below(anchors.len() as u64) as usize];
        // megablast, and two scorings that take the affine route
        let scoring = match case % 5 {
            0 => (2, -3, 5, 2),
            1 => (1, -3, 0, 0),
            _ => (1, -2, 0, 0),
        };
        let x_drop = 5 + rng.below(60) as i32;
        for traceback in [false, true] {
            let reference = run(
                &query,
                &subject,
                q_off,
                s_off,
                scoring,
                x_drop,
                &mut GreedyAlignScratch::new(),
                traceback,
            );
            let reused = run(
                &query, &subject, q_off, s_off, scoring, x_drop, &mut used, traceback,
            );
            assert_eq!(
                reused, reference,
                "case {case}: used scratch, traceback={traceback}"
            );
            if let Some(result) = &reference {
                aligned += 1;
                if let Some((ops, _)) = &result.7 {
                    if ops.iter().any(|op| *op != GapAlignOpType::Sub) {
                        with_gaps += 1;
                    }
                }
            }
        }
    }
    // the generator must produce real alignments, some of them gapped
    assert!(
        aligned > cases,
        "only {aligned} alignments in {cases} cases"
    );
    assert!(
        with_gaps * 20 > cases,
        "only {with_gaps} gapped alignments in {cases} cases"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:500-526,571-576
// if (fence_hit && *fence_hit) { return 0; }
// C writes into persistent scratch before the fence return.
#[test]
fn non_affine_fence_retains_written_scratch() {
    for traceback in [false, true] {
        for reverse in [false, true] {
            let mut mem = GreedyNonAffineMem::new();
            mem.ensure_base_capacity(14);
            mem.last_seq2_off[0].fill(123);
            let mut edit = GapPrelimEditBlock::new();
            let mut seed = GreedySeed::default();
            let mut fence = false;
            let (mut qlen, mut slen) = (0, 0);
            let score = blast_greedy_align(
                &[0, 1, 0],
                3,
                &[0, FENCE_SENTRY, 0],
                3,
                reverse,
                20,
                2,
                4,
                &mut qlen,
                &mut slen,
                4,
                &mut mem,
                traceback.then_some(&mut edit),
                &mut fence,
                &mut seed,
            );
            assert!(fence);
            assert_eq!((score, qlen, slen), (0, 1, 1));
            assert_eq!(mem.last_seq2_off[0][6], 1);
            for index in [4, 5, 7, 8] {
                assert_eq!(mem.last_seq2_off[0][index], INVALID_OFFSET);
            }
            assert_eq!(mem.last_seq2_off[0][3], 123);
        }
    }
}

#[test]
fn test_greedy_traceback_ncbi_regression_cases() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2762-2936
    let cases: &[(
        &[u8],
        &[u8],
        usize,
        usize,
        i32,
        i32,
        i32,
        i32,
        i32,
        (i32, i32, i32, i32, i32, usize),
        &[GapEditOp],
    )] = &[
        (
            &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            5,
            5,
            1,
            -2,
            0,
            0,
            54,
            (0, 11, 0, 11, 11, 11),
            &[GapEditOp::Sub(11)],
        ),
        (
            &[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            6,
            5,
            1,
            -2,
            0,
            0,
            54,
            (1, 12, 0, 11, 11, 11),
            &[GapEditOp::Sub(11)],
        ),
        (
            &[0, 1, 2, 3, 0, 1, 2],
            &[0, 1, 2, 0, 1, 2],
            3,
            3,
            1,
            -2,
            0,
            0,
            54,
            (0, 7, 0, 6, 3, 7),
            &[GapEditOp::Sub(3), GapEditOp::Ins(1), GapEditOp::Sub(3)],
        ),
        (
            &[0, 1, 2, 3, 0, 1],
            &[0, 1, 2, 3, 0, 1],
            3,
            3,
            1,
            -2,
            5,
            2,
            54,
            (0, 6, 0, 6, 6, 6),
            &[GapEditOp::Sub(6)],
        ),
        (
            &[0],
            &[0],
            0,
            0,
            1,
            -2,
            5,
            2,
            54,
            (0, 1, 0, 1, 1, 1),
            &[GapEditOp::Sub(1)],
        ),
    ];

    for (
        query,
        subject,
        q_off,
        s_off,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop,
        expected,
        expected_ops,
    ) in cases
    {
        let mut scratch = GreedyAlignScratch::new();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:479-490
        // Reuse the same scratch for the existing frozen cases and edits.
        for _ in 0..2 {
            let result = greedy_gapped_alignment_with_traceback(
                ContextQuery {
                    buf: query,
                    start: 0,
                    len: query.len() as i32,
                },
                AdjustedSubject {
                    buf: subject,
                    start: 0,
                    len: subject.len() as i32,
                },
                *q_off as i32,
                *s_off as i32,
                *reward,
                *penalty,
                *gap_open,
                *gap_extend,
                *x_drop,
                &mut scratch,
            )
            .expect("NCBI greedy traceback converges for the bounded fixture");
            assert_eq!(
                (result.0, result.1, result.2, result.3, result.4, result.5),
                *expected
            );
            assert_eq!(result.6, *expected_ops);
        }
    }
}

#[test]
fn test_greedy_reverse_and_empty_boundary_cases() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:380-453
    let reverse =
        greedy_align_one_direction_ex(&[0, 1, 2, 3], &[0, 1, 2, 3], 4, 4, 1, -2, 0, 0, 54, true);
    assert_eq!(reverse, (4, 4, 4, 4, 0, 0, 0));

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:471-472
    // Empty sequences: the early return, an alignment of length 0, in both directions.
    for reverse in [false, true] {
        assert_eq!(
            greedy_align_one_direction_ex(&[], &[], 0, 0, 1, -2, 0, 0, 54, reverse),
            (0, 0, 0, 0, 0, 0, 0)
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2783-2786,2839-2840
// A seed outside its context (remaining-gaps X1b-o, review 1 M2): past the context's end the
// right extension has a negative length and reads nothing, and the left one reads `q_off`
// letters back from the seed, through the sentinel (which matches nothing) into the context;
// before the context's start the left extension has a negative length, and the right one
// reads from the seed in the concatenated query up to the context's end. Each alignment here
// is one mismatch at the sentinel and 30 matches (score 30 - 2), whose query offsets lie
// outside the context as NCBI's do.
#[test]
fn a_greedy_seed_outside_its_context_aligns_on_the_concatenated_query() {
    let letters = |n: u8, a: u8| (0..n).map(|i| (i * 7 + a) % 4).collect::<Vec<u8>>();
    let mut scratch = GreedyAlignScratch::new();
    // Past the end: context 0 is 30 letters, the seed the second letter after it.
    let ctx0 = letters(30, 3);
    let ctx1 = letters(10, 0);
    let concat: Vec<u8> = [ctx0.as_slice(), &[15], ctx1.as_slice()].concat();
    let mut subject = ctx0.clone();
    subject.push((ctx0[29] + 2) % 4);
    subject.extend(ctx1.iter().map(|&b| (b + 1) % 4));
    let past_end = greedy_gapped_alignment_with_traceback(
        ContextQuery {
            buf: &concat,
            start: 0,
            len: 30,
        },
        AdjustedSubject {
            buf: &subject,
            start: 0,
            len: subject.len() as i32,
        },
        31,
        31,
        1,
        -2,
        0,
        0,
        54,
        &mut scratch,
    )
    .expect("NCBI aligns from a seed past its context's end");
    assert_eq!(
        (past_end.0, past_end.1, past_end.2, past_end.3, past_end.4),
        (0, 31, 0, 31, 28)
    );
    assert_eq!(past_end.6, vec![GapEditOp::Sub(31)]);
    // Before the start: context 1 is 30 letters after 10 of context 0 and the sentinel, the
    // seed on the sentinel (offset -1).
    let ctx0 = letters(10, 1);
    let ctx1 = letters(30, 2);
    let concat: Vec<u8> = [ctx0.as_slice(), &[15], ctx1.as_slice()].concat();
    let mut subject = vec![(ctx0[9] + 1) % 4, (ctx1[0] + 2) % 4];
    subject.extend_from_slice(&ctx1);
    let before_start = greedy_gapped_alignment_with_traceback(
        ContextQuery {
            buf: &concat,
            start: 11,
            len: 30,
        },
        AdjustedSubject {
            buf: &subject,
            start: 0,
            len: subject.len() as i32,
        },
        -1,
        1,
        1,
        -2,
        0,
        0,
        54,
        &mut scratch,
    )
    .expect("NCBI aligns from a seed before its context's start");
    assert_eq!(
        (
            before_start.0,
            before_start.1,
            before_start.2,
            before_start.3,
            before_start.4
        ),
        (-1, 30, 1, 32, 28)
    );
    assert_eq!(before_start.6, vec![GapEditOp::Sub(31)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2783-2787,2839-2840
// A start outside its subject range, which `AdjustSubjectRange` makes on a subject of 90000
// letters or more for a start more than 3000 letters outside its context (remaining-gaps
// X1b-o review 2 M1): with `s_avail` below 0 the right extension reads nothing and the left
// one reads `s_off` letters back from the start; with `s_off` below 0 the left extension
// reads nothing and the right one reads `s_avail` letters from the start, also when the
// range's pointer lies past the subject's end. The subject offsets are in the range.
#[test]
fn a_greedy_start_outside_its_subject_range_extends_to_the_other_side() {
    let query: Vec<u8> = (0..40usize).map(|i| ((i * 7 + 1) % 4) as u8).collect();
    let mut scratch = GreedyAlignScratch::new();
    let mut run = |subject: &[u8], start: usize, len: i32, q_off: i32, s_off: i32| {
        let r = greedy_gapped_alignment_with_traceback(
            ContextQuery {
                buf: &query,
                start: 0,
                len: 40,
            },
            AdjustedSubject {
                buf: subject,
                start,
                len,
            },
            q_off,
            s_off,
            1,
            -2,
            0,
            0,
            54,
            &mut scratch,
        )
        .expect("NCBI aligns from a start outside its subject range");
        assert_eq!(r.6.len(), 1);
        (r.0, r.1, r.2, r.3, r.4, r.6[0])
    };
    // `s_avail` = 10 - 20: five other letters, then the query; the range is the 10 letters
    // from the sixth, the start the query's 20th letter.
    let subject: Vec<u8> = [&[(query[0] + 1) % 4; 5][..], &query].concat();
    assert_eq!(
        run(&subject, 5, 10, 20, 20),
        (0, 20, 0, 20, 20, GapEditOp::Sub(20))
    );
    // The whole range from the same start aligns both ways.
    assert_eq!(
        run(&subject, 5, 40, 20, 20),
        (0, 40, 0, 40, 40, GapEditOp::Sub(40))
    );
    // `s_off` = -10: the range starts 30 letters into the query's copy, the start 10 letters
    // before it.
    assert_eq!(
        run(&query, 30, 10, 20, -10),
        (20, 40, -10, 10, 20, GapEditOp::Sub(20))
    );
    // The range's pointer 10 letters past the subject's end (length -10), the start 20
    // letters before it.
    assert_eq!(
        run(&query, 50, -10, 30, -20),
        (30, 40, -20, -10, 10, GapEditOp::Sub(10))
    );
}

/// A small xorshift generator for the random tests of this module.
struct TestRng(u64);

impl TestRng {
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

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:327-339,351-362
/// The affine kernel's scan (`affine_first_mismatch`) gives the byte loop's length and fence
/// flag (`find_first_mismatch_greedy`) on random pairs with ambiguous letters, fences and long
/// runs, in both directions, from starts inside, at and past the ends. Both slices run on past
/// their lengths by 0-9 bytes (bases that continue the other sequence, and fences), as the
/// reverse calls pass the whole query and subject with the offsets as lengths: a scan that read
/// past a length, or addressed a slice from its end instead of its length, would differ.
#[test]
fn affine_first_mismatch_equals_find_first_mismatch_greedy() {
    let mut rng = TestRng(0x0BAD_5EED_1234_5678);
    for _ in 0..20_000 {
        let len1 = rng.below(48) as usize;
        let seq1: Vec<u8> = (0..len1)
            .map(|_| {
                if rng.below(24) == 0 {
                    4 + rng.below(11) as u8
                } else {
                    rng.below(4) as u8
                }
            })
            .collect();
        let shift = rng.below(8) as usize;
        let len2 = rng.below(48) as usize;
        let seq2: Vec<u8> = (0..len2)
            .map(|j| match rng.below(40) {
                0 => FENCE_SENTRY,
                1..=4 => rng.below(4) as u8,
                _ => seq1.get(j + shift).copied().unwrap_or(0).min(3),
            })
            .collect();
        let s1 = rng.below(len1 as u64 + 3) as i32;
        let s2 = rng.below(len2 as u64 + 3) as i32;
        let (l1, l2) = (len1 as i32, len2 as i32);
        let (mut seq1, mut seq2) = (seq1, seq2);
        let (tail1, tail2) = (rng.below(10) as usize, rng.below(10) as usize);
        for i in 0..tail1 {
            // `seq2[j]` continues `seq1[j + shift]` (above)
            let b = (len1 + i)
                .checked_sub(shift)
                .and_then(|j| seq2.get(j))
                .copied()
                .unwrap_or(0)
                .min(3);
            seq1.push(if rng.below(3) == 0 {
                rng.below(4) as u8
            } else {
                b
            });
        }
        for i in 0..tail2 {
            let b = seq1.get(len2 + i + shift).copied().unwrap_or(0).min(3);
            seq2.push(match rng.below(12) {
                0 => FENCE_SENTRY,
                1..=3 => rng.below(4) as u8,
                _ => b,
            });
        }
        for reverse in [false, true] {
            let mut fence_old = false;
            let old =
                find_first_mismatch_greedy(&seq1, &seq2, l1, l2, s1, s2, reverse, &mut fence_old);
            let mut fence_new = false;
            let new = if reverse {
                affine_first_mismatch::<true>(&seq1, &seq2, l1, l2, s1, s2, &mut fence_new)
            } else {
                affine_first_mismatch::<false>(&seq1, &seq2, l1, l2, s1, s2, &mut fence_new)
            };
            assert_eq!(
                (new, fence_new),
                (old, fence_old),
                "reverse {reverse} s1 {s1} s2 {s2}"
            );
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:755-761,1000-1047
/// `reverse` reads both sequences from their ends, so the affine alignment of a pair with
/// `reverse` equals the forward alignment of the reversed pair: score, lengths, seed, fence
/// and edit script, on random related pairs with fences, ambiguous letters, several scoring
/// systems, X-drops from small to past every cell, short `max_dist` (no convergence), with and
/// without traceback. It runs both directions of the kernel (`blast_affine_greedy_align_dir`)
/// and its once-per-distance row checks (`prior_row_base`) on every distance.
#[test]
fn affine_greedy_reverse_equals_forward_on_reversed_pair() {
    type Out = (
        i32,
        i32,
        i32,
        bool,
        (i32, i32, i32),
        Vec<(GapAlignOpType, i32)>,
    );
    fn run(
        seq1: &[u8],
        seq2: &[u8],
        reverse: bool,
        scoring: (i32, i32, i32, i32),
        x_drop: i32,
        max_dist: i32,
        traceback: bool,
        mem: &mut (GreedyAffineMem, GreedyNonAffineMem),
    ) -> Out {
        let (mut len_1, mut len_2) = (0, 0);
        let mut fence_hit = false;
        let mut seed = GreedySeed::default();
        let mut block = GapPrelimEditBlock::new();
        let score = blast_affine_greedy_align(
            seq1,
            seq1.len() as i32,
            seq2,
            seq2.len() as i32,
            reverse,
            x_drop,
            scoring.0,
            scoring.1,
            scoring.2,
            scoring.3,
            &mut len_1,
            &mut len_2,
            max_dist,
            &mut mem.0,
            &mut mem.1,
            traceback.then_some(&mut block),
            &mut fence_hit,
            &mut seed,
        );
        let ops = block.edit_ops[..block.num_ops]
            .iter()
            .map(|op| (op.op_type, op.num))
            .collect();
        (
            score,
            len_1,
            len_2,
            fence_hit,
            (seed.start_q, seed.start_s, seed.match_length),
            ops,
        )
    }
    let mut rng = TestRng(0x1F2E_3D4C_5B6A_7988);
    let scorings = [(1, 2), (1, 3), (2, 3), (1, 1), (2, 5), (4, 5)];
    let x_drops = [2, 5, 10, 30, 100, 1_000, 20_000];
    let (mut fences, mut runs) = (0, 0);
    let mut mem = (GreedyAffineMem::new(), GreedyNonAffineMem::new());
    for _ in 0..600 {
        let len1 = 1 + rng.below(160) as usize;
        let seq1: Vec<u8> = (0..len1)
            .map(|_| {
                if rng.below(60) == 0 {
                    4 + rng.below(11) as u8
                } else {
                    rng.below(4) as u8
                }
            })
            .collect();
        let mut seq2 = Vec::new();
        for &b in &seq1 {
            match rng.below(100) {
                0..=7 => seq2.push(rng.below(4) as u8),
                8..=9 => {}
                10..=11 => {
                    seq2.push(b.min(3));
                    seq2.push(rng.below(4) as u8);
                }
                _ => seq2.push(b.min(3)),
            }
        }
        for _ in 0..rng.below(40) {
            seq2.push(rng.below(4) as u8);
        }
        if seq2.is_empty() {
            seq2.push(0);
        }
        if rng.below(4) == 0 {
            let at = rng.below(seq2.len() as u64) as usize;
            seq2[at] = FENCE_SENTRY;
        }
        let (reward, penalty) = scorings[rng.below(scorings.len() as u64) as usize];
        let gap_open = rng.below(7) as i32;
        let gap_extend = 1 + rng.below(4) as i32;
        let x_drop = x_drops[rng.below(x_drops.len() as u64) as usize];
        let max_dist = if rng.below(5) == 0 {
            1 + rng.below(8) as i32
        } else {
            (len1.max(seq2.len()) / 2 + 10) as i32
        };
        let scoring = (reward, penalty, gap_open, gap_extend);
        let rev1: Vec<u8> = seq1.iter().rev().copied().collect();
        let rev2: Vec<u8> = seq2.iter().rev().copied().collect();
        for traceback in [false, true] {
            let backward = run(
                &seq1, &seq2, true, scoring, x_drop, max_dist, traceback, &mut mem,
            );
            let forward = run(
                &rev1, &rev2, false, scoring, x_drop, max_dist, traceback, &mut mem,
            );
            assert_eq!(
                backward, forward,
                "scoring {scoring:?} x_drop {x_drop} max_dist {max_dist} traceback {traceback}"
            );
            fences += usize::from(forward.3);
            runs += 1;
        }
    }
    assert!(fences >= 20 && runs == 1_200, "fences {fences} runs {runs}");
}
