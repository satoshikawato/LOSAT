use super::x_kmerge::XUnitSort;
use super::*;
use crate::algorithm::tblastx::lookup::x_pair_lut::{
    x_assert_same_contexts, x_assert_same_group, x_assert_stub_scalars,
};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:343-346
// No computation: compares the direct groups with the groups `build_groups` filtered out of the
// reference builder's table (every field of every table, `q_min`, `q_max`), the stub with that
// table, and the contexts; the first difference stops the test.
fn lut_shadow_compare(
    lut: &XLutBatch,
    groups: &[XGroup],
    sliced: Option<&[XSliced]>,
    lookup: &BlastAaLookupTable,
    contexts: &[QueryContext],
) {
    assert_eq!(
        lut.groups.len(),
        groups.len(),
        "number of groups differs (direct vs reference)"
    );
    for (g, (direct, reference)) in lut.groups.iter().zip(groups).enumerate() {
        x_assert_same_group(
            g,
            &direct.lookup,
            (direct.q_min, direct.q_max),
            &reference.lookup,
            (reference.q_min, reference.q_max),
        );
    }
    // EXPERIMENT (LOSAT_X_TBXQSLICE): the slices, table by table, and each group's bounds.
    x_tbx_qslice::lut_shadow_compare(lut.sliced.as_deref(), sliced);
    let (stub, stub_contexts) = lut
        .stub
        .as_ref()
        .expect("the batch keeps the stub for the comparison");
    x_assert_stub_scalars(stub, lookup);
    x_assert_same_contexts(stub_contexts, contexts);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844
// No computation: compares the subject's combined list of the unit path with the reference
// loop's, every field, and stops at the first difference.
/// Check of one subject: the unit path's `combined_ungapped_hits` against the reference's.
fn shadow_compare(s_idx: usize, unit: &[UngappedHit], reference: &[UngappedHit]) {
    let n = unit.len().max(reference.len());
    for k in 0..n {
        match (unit.get(k), reference.get(k)) {
            (Some(a), Some(b)) if x_same_hit(a, b) => {}
            (a, b) => panic!(
                "subject {s_idx} frame {:?} index {k} differs (unit {} hits, reference {} hits): unit {a:?} vs reference {b:?}",
                a.or(b).map(|h| h.s_f_idx),
                unit.len(),
                reference.len()
            ),
        }
    }
}

/// A frame of amino acids (NCBISTDAA codes) with its sentinels.
fn frame(frame: i8, residues: &[u8]) -> QueryFrame {
    let mut aa_seq = Vec::with_capacity(residues.len() + 2);
    aa_seq.push(0);
    aa_seq.extend_from_slice(residues);
    aa_seq.push(0);
    QueryFrame {
        frame,
        aa_seq,
        aa_seq_nomask: None,
        aa_len: residues.len(),
        orig_len: residues.len() * 3 + 2,
        seg_masks: Vec::new(),
    }
}

/// NCBISTDAA codes of a letter string.
fn aa(s: &str) -> Vec<u8> {
    s.bytes()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect()
}

/// xorshift64*, for reproducible random residues.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// Random residues from a small alphabet (many repeated words and extensions).
fn random_residues(rng: &mut Rng, len: usize, alphabet: &[u8]) -> Vec<u8> {
    (0..len)
        .map(|_| alphabet[rng.below(alphabet.len() as u64) as usize])
        .collect()
}

fn lookup_for(frames: &[QueryFrame]) -> (BlastAaLookupTable, Vec<QueryContext>) {
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    build_ncbi_lookup(&[frames.to_vec()], 11, &params, true)
}

fn query_length_of(contexts: &[QueryContext]) -> i32 {
    contexts
        .last()
        .map(|c| c.frame_base + c.aa_len as i32)
        .unwrap_or(0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
// LOSAT_X_TBXCTXTABLE: the context the units' two-hit test takes (`x_context_of` on the
// batch's `x_ctx_table`) against the binary search `=0` runs (`get_context_idx`), for every
// offset from below the first frame base to past the batch's end, on 240 random layouts:
// 0 to 300 contexts, empty contexts (equal frame bases), a first frame base above 0, contexts
// of 1 to 9,000 letters and batches over 16,384 letters (coarse blocks of several offsets).
#[test]
fn the_context_table_gives_the_binary_search_index() {
    let frames = [frame(1, &aa("MKVLAAGIVRSTWYHE"))];
    let (mut lookup, _) = lookup_for(&frames);
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for layout in 0..240 {
        let contexts = rng.below(if layout % 3 == 0 { 300 } else { 12 }) as usize;
        let mut bases = Vec::with_capacity(contexts);
        let mut at = if layout % 5 == 0 {
            rng.below(50) as i32
        } else {
            0
        };
        for _ in 0..contexts {
            bases.push(at);
            at += match rng.below(6) {
                0 => 0,
                1 => 1 + rng.below(3) as i32,
                2 => 1 + rng.below(9_000) as i32,
                _ => 1 + rng.below(400) as i32,
            };
        }
        lookup.frame_bases = bases;
        lookup.num_contexts = contexts;
        lookup.query_length = at + rng.below(3) as i32;
        let table = x_ctx_table(&lookup, true);
        let off = x_ctx_table(&lookup, false);
        for q in -3..lookup.query_length + 3 {
            let reference = lookup.get_context_idx(q);
            if q >= 0 {
                let q = q as u32;
                assert_eq!(
                    x_context_of(&table, &lookup, q),
                    reference,
                    "layout {layout}, offset {q}: {:?}",
                    lookup.frame_bases
                );
                assert_eq!(x_context_of(&off, &lookup, q), reference);
            }
            // The table answers every offset of the batch itself; with `=0` none.
            if (0..lookup.query_length).contains(&q) {
                assert_eq!(table.get(q), Some(reference), "layout {layout}, offset {q}");
            }
            assert_eq!(off.get(q), None, "layout {layout}, offset {q}: `=0` table");
        }
    }
}

fn diag_array_size_for(query_length: i32, window: i32) -> i32 {
    let mut size: i32 = 1;
    while size < query_length + window {
        size <<= 1;
    }
    size
}

/// The reference two-hit stage of one subject on one masked table shared by the frames (the
/// frame loop, `frame_chunks::search_subject_frames`, with `scan_chunk_plain`): the batch
/// lookup, `(q - s) & mask`, the same per-hit body, the diagonal offset advanced after every
/// frame (`Blast_ExtendWordExit`). Returns the sorted init-HSP list of every frame, how often
/// each branch was taken, and how many hits met a cell whose flag 1 was written by a right
/// extension on another raw diagonal (the masked-aliasing case of module doc item 3). The
/// scores and the word length are those of `sc` (the extensions of the reference loop without
/// the sentinel walks: `extend_hit_two_hit_fitting` for BLOSUM62's instance,
/// `extend_hit_two_hit_with` for the other).
#[allow(clippy::too_many_arguments)]
fn masked_reference<M: TbxScore>(
    sc: M,
    lookup: &BlastAaLookupTable,
    contexts: &[QueryContext],
    s_frames: &[QueryFrame],
    cutoff_scores: &[i32],
    x_dropoff_per_context: &[i32],
    window: i32,
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
) -> (Vec<Vec<InitHSP>>, Vec<XHit>, usize) {
    let wordsize = sc.word() as i32;
    let mask = (diag_array.len() - 1) as u32;
    let mut per_frame = Vec::new();
    let mut outcomes = Vec::new();
    // The raw diagonal of the right extension that last set each cell's flag.
    let mut flag_writer = vec![i64::MIN; diag_array.len()];
    let mut aliased = 0usize;
    let array_size = OFFSET_ARRAY_SIZE + lookup.longest_chain.max(0);
    let mut offset_pairs = vec![OffsetPair::default(); array_size as usize];
    for (f, s_frame) in s_frames.iter().enumerate() {
        let subject = &s_frame.aa_seq[1..s_frame.aa_seq.len() - 1];
        let mut init_hsps = Vec::new();
        if subject.len() >= wordsize as usize {
            let failure = crate::algorithm::tblastx::subjects::XReadFailure::new();
            let ctx_of = x_ctx_table(lookup, true);
            let h = XHitCtx {
                diag_offset: *diag_offset,
                window,
                wordsize,
                lookup,
                ctx_of: &ctx_of,
                contexts,
                cutoff_scores,
                x_dropoff_per_context,
                subject,
                subject_sent: &[],
                s_f_idx: f,
                s_frame: s_frame.frame,
                s_idx: 0,
                s_len: 999,
                sc,
                query_batch: QueryBatch::new(contexts, &failure),
            };
            let seq_ranges = [(0, subject.len() as i32)];
            let mut scan_range = [0, 0, subject.len() as i32 - wordsize];
            while scan_range[1] <= scan_range[2] {
                let hits = s_blast_aa_scan_subject(
                    lookup,
                    subject,
                    &seq_ranges,
                    &mut offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                for pair in &offset_pairs[..hits as usize] {
                    let cell = (pair.q_off.wrapping_sub(pair.s_off) & mask) as usize;
                    let raw = pair.q_off as i64 - pair.s_off as i64;
                    if diag_array[cell].flag() != 0
                        && flag_writer[cell] != i64::MIN
                        && flag_writer[cell] != raw
                    {
                        aliased += 1;
                    }
                    let outcome = x_hit::<_, _, _, 0>(
                        &mut diag_array[cell],
                        pair.q_off,
                        pair.s_off,
                        &h,
                        &mut init_hsps,
                        &mut XNoProbe,
                    );
                    if matches!(outcome, XHit::Extended | XHit::Saved)
                        && diag_array[cell].flag() != 0
                    {
                        flag_writer[cell] = raw;
                    }
                    outcomes.push(outcome);
                }
            }
        }
        advance_tblastx_diag_offset(diag_offset, diag_array, window, subject.len());
        sort_init_hsps_by_score_ncbi(&mut init_hsps);
        per_frame.push(init_hsps);
    }
    (per_frame, outcomes, aliased)
}

/// The unit path's sorted init-HSP list of every frame (one chunk per frame here).
fn unit_lists<M: TbxScore>(pp: &XPairPar<M>, inp: &XSubject<'_>) -> Vec<Vec<InitHSP>> {
    let mut out = Vec::new();
    let mut scratch = XScratch::new(false);
    for (f, s_frame) in inp.s_frames_preliminary.iter().enumerate() {
        let chunk = SubjectChunk {
            offset: 0,
            length: s_frame.aa_len,
            overlap: 0,
        };
        let mut init_hsps = Vec::new();
        // The reference finds no hits in a frame shorter than a word.
        if chunk.length >= inp.wordsize as usize {
            for (g, group) in pp.groups.iter().enumerate() {
                if group.lookup.longest_chain > 0 {
                    let unit = pp.x_unit(inp, f, 0, chunk, g);
                    init_hsps.extend(pp.run_unit::<false>(&unit, inp, &mut scratch));
                }
            }
        }
        sort_init_hsps_by_score_ncbi(&mut init_hsps);
        out.push(init_hsps);
    }
    out
}

struct Case {
    frames: Vec<QueryFrame>,
    s_frames: Vec<QueryFrame>,
}

/// Compares the unit path with the masked reference started from `start` cells at
/// `start_offset`; returns the reference's branch outcomes.
fn check_case(case: &Case, start: DiagStruct, start_offset: i32, stale: bool) -> Vec<XHit> {
    let window = 40;
    let (lookup, contexts) = lookup_for(&case.frames);
    let query_length = query_length_of(&contexts);
    let size = diag_array_size_for(query_length, window);
    let cutoff_scores = vec![12; contexts.len()];
    let x_dropoff_per_context = vec![16; contexts.len()];
    let pp = XPairPar {
        groups: build_groups(&lookup, &contexts, X_MAX_GROUPS),
        stub: false,
        sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
        sc: Blosum62Score,
        sliced: None,
        mask_rule: None,
        extents: Vec::new(),
        ctx_of: x_ctx_table(&lookup, true),
    };
    let mut diag_array = vec![start; size as usize];
    let mut diag_offset = start_offset;
    if stale {
        // Leave the values of an earlier run of the same subject in the table.
        let _ = masked_reference(
            Blosum62Score,
            &lookup,
            &contexts,
            &case.s_frames,
            &cutoff_scores,
            &x_dropoff_per_context,
            window,
            &mut diag_array,
            &mut diag_offset,
        );
    }
    let (reference, outcomes, _) = masked_reference(
        Blosum62Score,
        &lookup,
        &contexts,
        &case.s_frames,
        &cutoff_scores,
        &x_dropoff_per_context,
        window,
        &mut diag_array,
        &mut diag_offset,
    );
    let inp = XSubject {
        s_frames_preliminary: &case.s_frames,
        s_frames: &case.s_frames,
        contexts: &contexts,
        cutoff_scores: &cutoff_scores,
        x_dropoff_per_context: &x_dropoff_per_context,
        window,
        wordsize: 3,
        s_idx: 0,
        s_len: 999,
        failure: &crate::algorithm::tblastx::subjects::XReadFailure::new(),
        passed: None,
        subject_besthit: false,
    };
    let units = unit_lists(&pp, &inp);
    assert_eq!(
        units, reference,
        "unit lists differ from the masked reference"
    );
    outcomes
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:595-599
// Line RG item X2c: frame 1 is R (7,000 letters, about 40,600 against itself) and frame 2 is
// P, the subject R, one letter, P, so R/R and P/P lie on one raw diagonal across the NULLB
// between the frames. With a raw X-drop of 40,000 the R/R extension's right walk goes on past
// the NULLB over P/P without scoring higher (about 40,600 - 32,768 + P/P), and its stop masks
// every P/P hit in the reference's one table. With the two frames in two groups, the unit of
// frame 2 alone extends P/P, whose left walk crosses the NULLB into R and scores higher there:
// a refusal, and a P/P HSP, that the reference does not make. `run_subject_units` runs the
// chunk's units again in group order on one table (`x_units_in_order`): the reference's list
// and no refusal, on 1 and 4 threads, with unmasked and masked unit tables and with frame 1's
// group cut into slices (LOSAT_X_TBXQSLICE), and the re-run alone into today's full lists.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn a_walk_past_a_nullb_into_another_group_is_seen_by_its_units() {
    let mut rng = Rng(0x5eed_02c3);
    let alphabet = aa("ACDEFGHIKLMNPQRSTVWY");
    let r = random_residues(&mut rng, 7_000, &alphabet);
    let p = random_residues(&mut rng, 200, &alphabet);
    let frames = vec![frame(1, &r), frame(2, &p)];
    let mut subject = r.clone();
    subject.push(alphabet[0]);
    subject.extend_from_slice(&p);
    let s_frames = vec![frame(1, &subject)];
    let window = 40;
    let (lookup, contexts) = lookup_for(&frames);
    let query_length = query_length_of(&contexts);
    let cutoff_scores = vec![30; contexts.len()];
    let x_dropoff_per_context = vec![40_000; contexts.len()];
    let mut diag_array =
        vec![DiagStruct::default(); diag_array_size_for(query_length, window) as usize];
    let mut diag_offset = window;
    let (reference, _, _) = masked_reference(
        Blosum62Score,
        &lookup,
        &contexts,
        &s_frames,
        &cutoff_scores,
        &x_dropoff_per_context,
        window,
        &mut diag_array,
        &mut diag_offset,
    );
    // The reference keeps R/R (frame 1, from the subject's start) and has no P/P HSP.
    let p_start = contexts[1].frame_base;
    assert!(reference[0]
        .iter()
        .any(|hsp| hsp.ctx_idx == 0 && hsp.s_start == 0 && hsp.score > 32_768));
    assert!(!reference[0]
        .iter()
        .any(|hsp| hsp.ctx_idx == 1 && hsp.q_start_absolute - p_start + 7_001 == hsp.s_start));
    let mut expect: Vec<UngappedHit> =
        get_ungapped_hsp_list(reference[0].clone(), &contexts, &s_frames);
    sort_ungapped_hits_by_score_ncbi(&mut expect);
    for (masked, sliced) in [(false, false), (true, false), (false, true)] {
        // Frame 1's group in seven slices of 1,000 offsets, frame 2's in one.
        let parts: Vec<Vec<i32>> = vec![(0..7).map(|i| i * 1_000).collect(), vec![p_start]];
        let pp = XPairPar {
            groups: if sliced {
                Vec::new()
            } else {
                build_groups(&lookup, &contexts, X_MAX_GROUPS)
            },
            stub: false,
            sent: false,
            sc: Blosum62Score,
            sliced: sliced.then(|| x_tbx_qslice::build_sliced(&lookup, &parts)),
            mask_rule: masked.then_some(XMaskRule {
                max_cells: usize::MAX,
                ratio: 1,
                band_cells: None,
            }),
            extents: if masked {
                x_group_extents(&contexts, X_MAX_GROUPS)
            } else {
                Vec::new()
            },
            ctx_of: x_ctx_table(&lookup, true),
        };
        assert_eq!(pp.x_group_count(), 2);
        assert!((0..2).all(|g| pp.x_group_has_entries(g)));
        let alone = crate::algorithm::tblastx::subjects::XReadFailure::new();
        let failure = crate::algorithm::tblastx::subjects::XReadFailure::new();
        let inp = XSubject {
            s_frames_preliminary: &s_frames,
            s_frames: &s_frames,
            contexts: &contexts,
            cutoff_scores: &cutoff_scores,
            x_dropoff_per_context: &x_dropoff_per_context,
            window,
            wordsize: 3,
            s_idx: 0,
            // As the reference's (`masked_reference`).
            s_len: 999,
            failure: &alone,
            passed: None,
            subject_besthit: false,
        };
        let chunk = SubjectChunk {
            offset: 0,
            length: subject.len(),
            overlap: 0,
        };
        let units: Vec<XUnit> = (0..2).map(|g| pp.x_unit(&inp, 0, 0, chunk, g)).collect();
        // Each unit alone: the refusal and the P/P HSP the reference does not make.
        if !sliced {
            assert_ne!(unit_lists(&pp, &inp), reference, "masked {masked}");
            assert!(
                alone.take().is_some(),
                "masked {masked}: the units alone refuse"
            );
        }
        let inp = XSubject {
            failure: &failure,
            ..inp
        };
        // The re-run alone, into today's lists: the reference's sorted list.
        let mut again = XUnitOutputs::Full(vec![Vec::new(), Vec::new()]);
        pp.x_units_in_order(&units, &inp, &mut again);
        let XUnitOutputs::Full(lists) = again else {
            unreachable!("the form given")
        };
        let mut joined: Vec<InitHSP> = lists.into_iter().flatten().collect();
        sort_init_hsps_by_score_ncbi(&mut joined);
        assert_eq!(
            joined, reference[0],
            "masked {masked} sliced {sliced}: the re-run"
        );
        assert!(
            failure.take().is_none(),
            "sliced {sliced}: the re-run refuses nothing"
        );
        for threads in [1, 4] {
            let merged = crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
                Ok(XPairPar::<Blosum62Score>::merge_subject(
                    pp.run_subject_units(&inp, pool),
                    &inp,
                ))
            })
            .expect("search pool");
            assert!(
                failure.take().is_none(),
                "masked {masked} sliced {sliced}, {threads} threads: no refusal"
            );
            assert_eq!(
                merged.len(),
                expect.len(),
                "masked {masked} sliced {sliced}, {threads} threads"
            );
            for (got, want) in merged.iter().zip(&expect) {
                assert!(
                    x_same_hit(got, want),
                    "masked {masked} sliced {sliced}, {threads} threads"
                );
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:337-360
#[test]
fn sub_lookup_chains_remerge_to_the_full_chain() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let alphabet = aa("ACDEKLMNPQRSTVWY");
    let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
        .iter()
        .map(|&f| frame(f, &random_residues(&mut rng, 300, &alphabet)))
        .collect();
    let (lookup, contexts) = lookup_for(&frames);
    let groups = build_groups(&lookup, &contexts, X_MAX_GROUPS);
    assert_eq!(groups.len(), 6);
    // The tables of 300-letter contexts are sparse: with LOSAT_X_TBXSCANBF they keep only
    // their cells with entries (`XGroup::new`), read here through `XGroup::hits`.
    assert_eq!(
        groups.iter().all(|g| !g.rank.is_empty()),
        crate::algorithm::tblastx::x_scan_bf::x_tbx_scan_bf()
    );
    for idx in 0..lookup.backbone.len() {
        let full = lookup.get_hits(idx);
        let rank = |q: i32| {
            full.iter()
                .position(|&x| x == q)
                .expect("entry of the chain")
        };
        let mut merged: Vec<i32> = Vec::new();
        for (g, group) in groups.iter().enumerate() {
            let sub = group.hits(idx);
            assert_eq!(
                crate::algorithm::tblastx::lookup::pv_test(&group.lookup.pv, idx),
                !sub.is_empty()
            );
            assert!(sub.len() as i32 <= group.lookup.longest_chain);
            // Order kept within the sub-chain, all entries of context g.
            assert!(sub.windows(2).all(|w| rank(w[0]) < rank(w[1])));
            assert!(sub.iter().all(|&q| lookup.get_context_idx(q) == g
                && q as u32 >= group.q_min
                && q as u32 <= group.q_max));
            merged.extend_from_slice(sub);
        }
        merged.sort_by_key(|&q| rank(q));
        assert_eq!(merged, full, "cell {idx}");
    }
    let total: usize = groups
        .iter()
        .map(|g| {
            g.lookup
                .backbone
                .iter()
                .map(|c| c.num_used as usize)
                .sum::<usize>()
        })
        .sum();
    let full_total: usize = lookup.backbone.iter().map(|c| c.num_used as usize).sum();
    assert_eq!(total, full_total);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-61
#[test]
fn last_hit_fits_edges() {
    // Default window, a query of 20 M letters, the longest chunk NCBI makes: far inside.
    let ql = 20_000_000;
    let size = diag_array_size_for(ql, 40);
    assert!(last_hit_fits(5_000_000 / 3, ql, size));
    // A window just under 2^29 (the table has 2^30 cells): the reference's second scan
    // stores values of 2^30 and more, so no unit path.
    let window = 536_870_000;
    let size = diag_array_size_for(2_000, window);
    assert_eq!(size, 1 << 30);
    assert!(!last_hit_fits(20_000, 2_000, size));
    // The bound itself: slack + 2 * chunk + margin one below and at 2^29.
    let size = 1 << 29;
    let ql = (1 << 29) - (1 << 28);
    let half_margin = (X_LAST_HIT_MARGIN / 2) as usize;
    assert!(last_hit_fits((1 << 27) - half_margin - 1, ql, size));
    assert!(!last_hit_fits((1 << 27) - half_margin, ql, size));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-61
#[test]
fn guard_edges() {
    let window = 40;
    // Query length just under a power of two: slack = window, contexts far longer, long
    // subject: no unit path.
    let ql = (1 << 20) - window;
    let size = diag_array_size_for(ql, window);
    assert_eq!(size, 1 << 20);
    assert!(!guard(ql as usize / 6, 1_000_000, ql, size));
    // Same query, subject chunk within the slack: no aliasing, units exact.
    assert!(guard(ql as usize / 6, window as usize, ql, size));
    assert!(!guard(ql as usize / 6, window as usize + 1, ql, size));
    // One letter longer: the table doubles, the slack covers every context.
    let size2 = diag_array_size_for(ql + 1, window);
    assert_eq!(size2, 1 << 21);
    assert!(guard((ql as usize + 1) / 6, 5_000_000, ql + 1, size2));
    // Context exactly the slack long passes, one more letter fails.
    let slack = (size2 - (ql + 1)) as usize;
    assert!(guard(slack, 5_000_000, ql + 1, size2));
    assert!(!guard(slack + 1, 5_000_000, ql + 1, size2));
    // p03 and d04 (map §5).
    assert!(guard(95_687, 98_048, 574_123, 1 << 20));
    assert!(guard(220_702, 202_064, 1_324_217, 1 << 21));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560-573
#[test]
fn context_boundary_diagonal() {
    // Context 0 ends with WWW, context 1 starts with CCC: on the raw diagonal of the subject
    // "WWWGCCC..." the CCC hit of context 1 follows the WWW hit of context 0 by 4 letters
    // and fails the context check in the reference; in the unit it is the first hit of a
    // fresh cell. A second CCC word later on the same diagonal then extends in both.
    let mut c0 = aa("KLMNPQRSTVWW");
    c0.extend(aa("WWW"));
    let c1 = aa("CCCHHHIIIEEECCC");
    let frames = vec![frame(1, &c0), frame(2, &c1)];
    let s = aa("AAWWWWCCCHHHIIIEEECCCAA");
    let case = Case {
        frames,
        s_frames: vec![frame(1, &s)],
    };
    let outcomes = check_case(&case, DiagStruct::default(), 40, false);
    assert!(outcomes.contains(&XHit::ContextBoundary), "{outcomes:?}");
    assert!(outcomes.contains(&XHit::Saved), "{outcomes:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-104
#[test]
fn cleared_cells_hold_minus_window() {
    let window = 40;
    let cell = DiagStruct::clear(window);
    // Sign extension of the 31-bit field (round-1 bug: 2^31 - window).
    assert_eq!(cell.last_hit(), -window);
    assert_eq!(cell.flag(), 0);
    // A hit on a cleared cell at diag_offset = window is "too far", like a fresh unit cell.
    let mut rng = Rng(7);
    let alphabet = aa("ACDEKW");
    let q = random_residues(&mut rng, 60, &alphabet);
    let s = random_residues(&mut rng, 80, &alphabet);
    let case = Case {
        frames: vec![frame(1, &q), frame(-1, &q[10..])],
        s_frames: vec![frame(1, &s), frame(-1, &s[5..])],
    };
    let outcomes = check_case(&case, DiagStruct::clear(window), window, false);
    assert!(outcomes.contains(&XHit::Saved));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
#[test]
fn diag_offset_near_the_reset_limit() {
    let mut rng = Rng(11);
    let alphabet = aa("ACDEKW");
    let q = random_residues(&mut rng, 70, &alphabet);
    let s = random_residues(&mut rng, 90, &alphabet);
    let s_frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
        .iter()
        .map(|&f| frame(f, &s[(f.unsigned_abs() as usize)..]))
        .collect();
    let case = Case {
        frames: vec![frame(1, &q), frame(2, &q[3..]), frame(-1, &q[7..])],
        s_frames,
    };
    // Stale values below the limit, then frames that cross it (reset to window mid-subject).
    for start in [i32::MAX / 4 - 400, i32::MAX / 4 - 100, i32::MAX / 4] {
        check_case(&case, DiagStruct::default(), start, true);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
#[test]
fn random_cases_with_and_without_aliasing() {
    let mut rng = Rng(0x1234_5678_9abc_def1);
    let alphabets = [aa("ACDEKW"), aa("AW"), aa("ACDEFGHIKLMNPQRSTVWY")];
    let mut aliased = 0usize;
    let mut checked = 0usize;
    for case_no in 0..300 {
        let alphabet = &alphabets[case_no % alphabets.len()];
        let n_ctx = 1 + rng.below(6) as usize;
        let frames: Vec<QueryFrame> = (0..n_ctx)
            .map(|c| {
                let len = 3 + rng.below(40) as usize;
                frame(
                    [1i8, 2, 3, -1, -2, -3][c],
                    &random_residues(&mut rng, len, alphabet),
                )
            })
            .collect();
        let s_len = 3 + rng.below(400) as usize;
        let s_frames: Vec<QueryFrame> = (0..1 + rng.below(3) as usize)
            .map(|f| {
                frame(
                    [1i8, 2, 3][f],
                    &random_residues(&mut rng, s_len - f, alphabet),
                )
            })
            .collect();
        let (_, contexts) = lookup_for(&frames);
        let window = 40;
        let query_length = query_length_of(&contexts);
        let size = diag_array_size_for(query_length, window);
        let max_ctx = contexts.iter().map(|c| c.aa_len).max().unwrap_or(0);
        if !guard(max_ctx, s_len, query_length, size) {
            continue;
        }
        if s_len as i64 > (size - query_length) as i64 {
            aliased += 1;
        }
        let case = Case { frames, s_frames };
        let start = match case_no % 3 {
            0 => (DiagStruct::default(), window),
            1 => (DiagStruct::clear(window), window),
            _ => (DiagStruct::default(), i32::MAX / 4 - 50),
        };
        check_case(&case, start.0, start.1, case_no % 2 == 1);
        checked += 1;
    }
    assert!(checked > 100, "checked {checked}");
    assert!(aliased > 20, "aliased {aliased}");
}

// No NCBI counterpart: LOSAT_X_TBXBAND (module doc item 3). Batches the first guard rejects: a
// context and the subject longer than the slack, the subject made of stretches copied from the
// long context so that extensions run past `slack + 3 + (wordsize - 1)` letters, the bound the
// earlier argument could not exclude. The units' lists are the masked reference's in every
// case, also with the batched extensions of LOSAT_X_TBXEXTBATCH.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-528,598-604
#[test]
fn units_equal_the_masked_reference_where_the_first_guard_fails() {
    let mut rng = Rng(0x7e57_ab1e_5eed_0001);
    let alphabets = [aa("AW"), aa("ACW"), aa("ACDEKW")];
    let window = 40;
    let (mut checked, mut long_reach, mut aliased_cases) = (0usize, 0usize, 0usize);
    for case_no in 0..400 {
        let alphabet = &alphabets[case_no % alphabets.len()];
        // One long context (longer than the slack) and up to two short ones.
        let long = 90 + rng.below(120) as usize;
        let mut frames = vec![frame(1, &random_residues(&mut rng, long, alphabet))];
        for f in 0..rng.below(3) as usize {
            let len = 3 + rng.below(20) as usize;
            frames.push(frame(
                [2i8, 3][f],
                &random_residues(&mut rng, len, alphabet),
            ));
        }
        let s_len = 150 + rng.below(350) as usize;
        let mut residues = Vec::with_capacity(s_len);
        while residues.len() < s_len {
            if rng.below(2) == 0 {
                let from = rng.below(long as u64 / 2) as usize;
                let take = (long - from).min(s_len - residues.len());
                residues.extend_from_slice(&frames[0].aa_seq[1 + from..1 + from + take]);
            } else {
                let take = (5 + rng.below(60) as usize).min(s_len - residues.len());
                residues.extend(random_residues(&mut rng, take, alphabet));
            }
        }
        let s_frames = vec![frame(1, &residues)];
        let (lookup, contexts) = lookup_for(&frames);
        let query_length = query_length_of(&contexts);
        let size = diag_array_size_for(query_length, window);
        let max_ctx = contexts.iter().map(|c| c.aa_len).max().unwrap_or(0);
        if guard(max_ctx, s_len, query_length, size) {
            continue;
        }
        assert!(last_hit_fits(s_len, query_length, size));
        let slack = size - query_length;
        let cutoff_scores = vec![12; contexts.len()];
        let x_dropoff_per_context = vec![16; contexts.len()];
        let mut diag_array = vec![DiagStruct::default(); size as usize];
        let mut diag_offset = window;
        let (reference, outcomes, aliased) = masked_reference(
            Blosum62Score,
            &lookup,
            &contexts,
            &s_frames,
            &cutoff_scores,
            &x_dropoff_per_context,
            window,
            &mut diag_array,
            &mut diag_offset,
        );
        let inp = XSubject {
            s_frames_preliminary: &s_frames,
            s_frames: &s_frames,
            contexts: &contexts,
            cutoff_scores: &cutoff_scores,
            x_dropoff_per_context: &x_dropoff_per_context,
            window,
            wordsize: 3,
            s_idx: 0,
            s_len: 999,
            failure: &crate::algorithm::tblastx::subjects::XReadFailure::new(),
            passed: None,
            subject_besthit: false,
        };
        let pp = XPairPar {
            groups: build_groups(&lookup, &contexts, X_MAX_GROUPS),
            stub: false,
            sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
            sc: Blosum62Score,
            sliced: None,
            mask_rule: None,
            extents: Vec::new(),
            ctx_of: x_ctx_table(&lookup, true),
        };
        assert_eq!(
            unit_lists(&pp, &inp),
            reference,
            "case {case_no}: unit lists differ from the masked reference"
        );
        assert_eq!(
            unit_lists_batched(&pp, &inp),
            reference,
            "case {case_no}: batched unit lists differ from the masked reference"
        );
        // A saved extension that ran further right of its seed than the earlier argument's
        // bound allowed: the reach `s_last_off - s0` is at least `s_end - s_seed - 1`
        // (`s_end` is one past the best end, which is at most the last letter scanned), and the
        // earlier argument fails only for a reach above `slack + 3 + 2`.
        if reference
            .iter()
            .flatten()
            .any(|h| h.s_end - h.s_seed - 1 > slack + 3 + 2)
        {
            long_reach += 1;
        }
        // A hit of the reference met a flag set by a right extension on another raw diagonal
        // sharing its cell (the masked-aliasing case), and took the unmasked branch.
        if aliased > 0 {
            aliased_cases += 1;
        }
        assert!(outcomes.len() >= aliased);
        checked += 1;
    }
    assert!(checked > 150, "checked {checked}");
    assert!(long_reach > 50, "long_reach {long_reach}");
    assert!(aliased_cases > 20, "aliased_cases {aliased_cases}");
}

// No NCBI counterpart: LOSAT_X_TBXUNITMASK (module doc). Units on masked tables against the
// same units on unmasked tables, in every hit loop (with and without LOSAT_X_TBXPREFETCH's hint
// and LOSAT_X_TBXHITBF, and LOSAT_X_TBXEXTBATCH's batched extensions), and the sorted lists
// against the reference loop on its own masked table. Query groups far shorter than the
// subject, so that many raw diagonals share each masked cell; subjects made of pieces of the
// query, so that long right extensions set flags that hits of aliasing diagonals then meet.
#[test]
fn masked_unit_tables_give_the_unmasked_units() {
    let mut rng = Rng(0x7e57_ab1e_5eed_0b1b);
    let alphabets = [aa("AW"), aa("ACW"), aa("ACDEKW")];
    let frame_numbers = [1i8, 2, 3, -1, -2, -3];
    let (mut masked_units, mut aliased_cases, mut saved) = (0usize, 0usize, 0usize);
    let (mut exact_fits, mut wide_groups, mut control_differs) = (0usize, 0usize, 0usize);
    let mut sent_units = 0usize;
    // LOSAT_X_TBXUNGAPSENT is on (the batches' contexts have their sentinels).
    let on = x_tbx_sent_contexts(std::iter::empty());
    for case_no in 0..300 {
        let alphabet = &alphabets[case_no % alphabets.len()];
        let window = [40, 3 + x_tbx_extbatch::BATCH as i32, 64][case_no % 3];
        // Every fourth case two queries (7 to 12 contexts: groups of several contexts); every
        // fifth one context lengths whose group extent plus the window is a power of two (a
        // group's extent is its length and the sentinel after it, the last group's its length).
        let queries = if case_no % 4 == 3 { 2 } else { 1 };
        let n_frames = if queries == 2 {
            7 + rng.below(6) as usize
        } else {
            1 + rng.below(6) as usize
        };
        let exact = queries == 1 && case_no % 5 == 0;
        let mut frames = Vec::with_capacity(n_frames);
        for i in 0..n_frames {
            let len = if exact {
                let cells = [64usize, 128][i % 2].max((window as usize + 5).next_power_of_two());
                cells - window as usize - usize::from(i + 1 < n_frames)
            } else {
                4 + rng.below(80) as usize
            };
            frames.push(frame(
                frame_numbers[i % 6],
                &random_residues(&mut rng, len, alphabet),
            ));
        }
        let s_len = 1500 + rng.below(2500) as usize;
        let mut residues = Vec::with_capacity(s_len);
        while residues.len() < s_len {
            if rng.below(3) != 0 {
                let q = &frames[rng.below(n_frames as u64) as usize];
                let len = q.aa_seq.len() - 2;
                let from = rng.below(len as u64) as usize;
                let take = (len - from).min(s_len - residues.len());
                residues.extend_from_slice(&q.aa_seq[1 + from..1 + from + take]);
            } else {
                let take = (3 + rng.below(40) as usize).min(s_len - residues.len());
                residues.extend(random_residues(&mut rng, take, alphabet));
            }
        }
        let s_frames = vec![frame(1, &residues)];
        let (lookup, contexts) = if queries == 2 {
            let (a, b) = frames.split_at(n_frames / 2);
            let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
            build_ncbi_lookup(&[a.to_vec(), b.to_vec()], 11, &params, true)
        } else {
            lookup_for(&frames)
        };
        let query_length = query_length_of(&contexts);
        let size = diag_array_size_for(query_length, window);
        let cutoff_scores = vec![12; contexts.len()];
        let x_dropoff_per_context = vec![16; contexts.len()];
        let mut diag_array = vec![DiagStruct::default(); size as usize];
        let mut diag_offset = window;
        let (reference, _, _) = masked_reference(
            Blosum62Score,
            &lookup,
            &contexts,
            &s_frames,
            &cutoff_scores,
            &x_dropoff_per_context,
            window,
            &mut diag_array,
            &mut diag_offset,
        );
        let inp = XSubject {
            s_frames_preliminary: &s_frames,
            s_frames: &s_frames,
            contexts: &contexts,
            cutoff_scores: &cutoff_scores,
            x_dropoff_per_context: &x_dropoff_per_context,
            window,
            wordsize: 3,
            s_idx: 0,
            s_len: 999,
            failure: &crate::algorithm::tblastx::subjects::XReadFailure::new(),
            passed: None,
            subject_besthit: false,
        };
        let plain = XPairPar {
            groups: build_groups(&lookup, &contexts, X_MAX_GROUPS),
            stub: false,
            sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
            sc: Blosum62Score,
            sliced: None,
            mask_rule: None,
            extents: Vec::new(),
            ctx_of: x_ctx_table(&lookup, true),
        };
        // Ratio 1: every unit whose unmasked table is at least the masked one's size.
        let pp = XPairPar {
            groups: build_groups(&lookup, &contexts, X_MAX_GROUPS),
            stub: false,
            sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
            sc: Blosum62Score,
            sliced: None,
            mask_rule: Some(XMaskRule {
                max_cells: usize::MAX,
                ratio: 1,
                band_cells: None,
            }),
            extents: x_group_extents(&contexts, X_MAX_GROUPS),
            ctx_of: x_ctx_table(&lookup, true),
        };
        assert_eq!(pp.extents.len(), pp.groups.len());
        // Negative control: tables of half the masked size (below extent + window, against
        // the module doc's condition), where the units may differ.
        let control = XPairPar {
            groups: build_groups(&lookup, &contexts, X_MAX_GROUPS),
            stub: false,
            sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
            sc: Blosum62Score,
            sliced: None,
            mask_rule: pp.mask_rule,
            extents: pp
                .extents
                .iter()
                .map(|&e| (x_masked_cells(e, window) / 2).saturating_sub(window as usize))
                .collect(),
            ctx_of: x_ctx_table(&lookup, true),
        };
        let chunk = SubjectChunk {
            offset: 0,
            length: s_len,
            overlap: 0,
        };
        let mut scratch = XScratch::new(false);
        let mut aliased = false;
        for g in 0..pp.groups.len() {
            if pp.groups[g].lookup.longest_chain <= 0 {
                continue;
            }
            let unit = pp.x_unit(&inp, 0, 0, chunk, g);
            let mask = pp.x_unit_mask(&unit, window).expect("a masked table");
            assert!(plain.x_unit_mask(&unit, window).is_none());
            masked_units += 1;
            // The masked table covers the group's extent plus the window (module doc).
            assert!(mask as usize + 1 >= pp.extents[g] + window as usize);
            exact_fits += usize::from(mask as usize + 1 == pp.extents[g] + window as usize);
            wide_groups +=
                usize::from(group_starts(&contexts, X_MAX_GROUPS).len() < contexts.len());
            let mut want: Vec<InitHSP> = Vec::new();
            plain.scan_unit_hits::<_, XNoProbe, false, false, 0>(
                &unit,
                &inp,
                &mut scratch,
                &mut want,
                &mut XNoProbe,
            );
            saved += want.len();
            // Two saved extensions of the unit on raw diagonals that share a masked cell.
            let diagonals: Vec<i64> = want
                .iter()
                .map(|h| h.q_seed_absolute as i64 - h.s_seed as i64)
                .collect();
            let cells = mask as i64 + 1;
            aliased |= diagonals.iter().any(|&a| {
                diagonals
                    .iter()
                    .any(|&b| a != b && (a - b).rem_euclid(cells) == 0)
            });
            // LOSAT_X_TBXUNGAPSENT: the subject is one whole frame with its sentinels, so its
            // units take the sentinel walks whenever the batch does.
            let sent = pp.x_unit_sent(&unit, &inp);
            assert_eq!(sent, pp.sent);
            sent_units += usize::from(sent);
            let mut lists: Vec<(&str, Vec<InitHSP>)> = Vec::new();
            for name in [
                "scan",
                "scan pf",
                "bf",
                "bf pf",
                "scan sent",
                "bf sent",
                "bf pf sent",
            ] {
                if name.ends_with("sent") && !sent {
                    continue;
                }
                let mut got: Vec<InitHSP> = Vec::new();
                match name {
                    "scan" => pp.scan_unit_hits::<_, XNoProbe, false, false, 0>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        &mut XNoProbe,
                    ),
                    "scan pf" => pp.scan_unit_hits::<_, XNoProbe, true, false, 0>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        &mut XNoProbe,
                    ),
                    "bf" => pp.scan_unit_hits::<_, XNoProbe, false, true, 0>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        &mut XNoProbe,
                    ),
                    "bf pf" => pp.scan_unit_hits::<_, XNoProbe, true, true, 0>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        &mut XNoProbe,
                    ),
                    "scan sent" => pp.scan_unit_hits::<_, XNoProbe, false, false, 3>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        &mut XNoProbe,
                    ),
                    "bf sent" => pp.scan_unit_hits::<_, XNoProbe, false, true, 3>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        &mut XNoProbe,
                    ),
                    _ => pp.scan_unit_hits::<_, XNoProbe, true, true, 3>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        &mut XNoProbe,
                    ),
                }
                assert_eq!(scratch.diag.len(), mask as usize + 1);
                lists.push((name, got));
            }
            for (name, len) in [("batched", x_tbx_extbatch::BATCH), ("batched 1", 1)] {
                let mut got: Vec<InitHSP> = Vec::new();
                pp.scan_unit_batched_hits::<_, false, false, 0>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut got,
                    len,
                );
                lists.push((name, got));
                let mut got: Vec<InitHSP> = Vec::new();
                pp.scan_unit_batched_hits::<_, true, true, 0>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut got,
                    len,
                );
                lists.push((name, got));
                if sent {
                    let mut got: Vec<InitHSP> = Vec::new();
                    pp.scan_unit_batched_hits::<_, true, true, 3>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        len,
                    );
                    lists.push((name, got));
                    let mut got: Vec<InitHSP> = Vec::new();
                    pp.scan_unit_batched_hits::<_, false, false, 3>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        len,
                    );
                    lists.push((name, got));
                }
            }
            for (name, got) in &lists {
                assert_eq!(
                    got, &want,
                    "case {case_no} group {g}: {name} on the masked table differs"
                );
            }
            if control.x_unit_mask(&unit, window).map(|m| m as usize + 1)
                == Some((mask as usize).div_ceil(2))
            {
                let mut got: Vec<InitHSP> = Vec::new();
                control.scan_unit_hits::<_, XNoProbe, false, false, 0>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut got,
                    &mut XNoProbe,
                );
                control_differs += usize::from(got != want);
            }
        }
        assert_eq!(
            unit_lists(&pp, &inp),
            reference,
            "case {case_no}: masked unit lists differ from the reference"
        );
        assert_eq!(
            unit_lists_batched(&pp, &inp),
            reference,
            "case {case_no}: masked batched unit lists differ from the reference"
        );
        aliased_cases += usize::from(aliased);
    }
    assert!(masked_units > 600, "masked_units {masked_units}");
    assert_eq!(sent_units, if on { masked_units } else { 0 });
    assert!(aliased_cases > 100, "aliased_cases {aliased_cases}");
    assert!(saved > 10_000, "saved {saved}");
    assert!(exact_fits > 50, "exact_fits {exact_fits}");
    assert!(wide_groups > 100, "wide_groups {wide_groups}");
    assert!(control_differs > 10, "control_differs {control_differs}");
}

// No NCBI counterpart: the sizes of LOSAT_X_TBXUNITMASK (module doc).
#[test]
fn masked_table_sizes_and_rule() {
    assert_eq!(x_masked_cells(0, 40), 64);
    assert_eq!(x_masked_cells(24, 40), 64);
    assert_eq!(x_masked_cells(25, 40), 128);
    assert_eq!(x_masked_cells(1 << 20, 0), 1 << 20);
    assert_eq!(x_masked_cells((1 << 20) - 40, 40), 1 << 20);
    assert_eq!(x_masked_cells((1 << 20) - 39, 40), 1 << 21);
    // Small masked tables where the unmasked one has at least `ratio` times their cells; in a
    // band batch, every unit whose chunk is longer than the reference table.
    let rule = |max_cells: usize, band_cells: Option<i32>| XMaskRule {
        max_cells,
        ratio: 4,
        band_cells,
    };
    assert_eq!(x_unit_mask_for(&rule(64, None), 24, 40, 256, 200), Some(63));
    assert_eq!(x_unit_mask_for(&rule(64, None), 24, 40, 255, 200), None);
    assert_eq!(x_unit_mask_for(&rule(63, None), 24, 40, 256, 200), None);
    assert_eq!(
        x_unit_mask_for(&rule(128, None), 25, 40, 512, 200),
        Some(127)
    );
    assert_eq!(x_unit_mask_for(&rule(128, None), 25, 40, 511, 200), None);
    assert_eq!(x_unit_mask_for(&rule(0, Some(256)), 24, 40, 300, 256), None);
    assert_eq!(
        x_unit_mask_for(&rule(0, Some(256)), 24, 40, 300, 257),
        Some(63)
    );
    assert_eq!(
        x_unit_mask_for(&rule(0, None), 24, 40, 1 << 20, 1 << 20),
        None
    );
    // The size limit by thread count: 8,192 cells at one, 4,096 at 2-7, the switch's KiB at 8+.
    assert_eq!(x_mask_max_cells(256, 1), 8192);
    assert_eq!(x_mask_max_cells(2, 1), 512);
    assert_eq!(x_mask_max_cells(256, 2), 4096);
    assert_eq!(x_mask_max_cells(256, 7), 4096);
    assert_eq!(x_mask_max_cells(256, 8), 65536);
    assert_eq!(x_mask_max_cells(2, 4), 512);
    assert_eq!(x_mask_max_cells(2, 16), 512);
    // The groups: one per context up to X_MAX_GROUPS, from a context's first offset to the
    // next group's (sentinels included), the last to the end of its last context.
    let frames: Vec<QueryFrame> = [5usize, 7, 9]
        .iter()
        .zip([1i8, 2, 3])
        .map(|(&len, f)| frame(f, &vec![aa("A")[0]; len]))
        .collect();
    let (_, contexts) = lookup_for(&frames);
    let extents = x_group_extents(&contexts, X_MAX_GROUPS);
    assert_eq!(extents.len(), 3);
    assert_eq!(
        extents[0],
        (contexts[1].frame_base - contexts[0].frame_base) as usize
    );
    assert_eq!(
        extents[1],
        (contexts[2].frame_base - contexts[1].frame_base) as usize
    );
    assert_eq!(extents[2], 9);
    assert!(extents[0] > 5 && extents[1] > 7);
}

/// The unit lists of `inp` with the batched extensions of LOSAT_X_TBXEXTBATCH (`unit_lists`
/// with `run_unit_batched`).
fn unit_lists_batched<M: TbxScore>(pp: &XPairPar<M>, inp: &XSubject<'_>) -> Vec<Vec<InitHSP>> {
    let mut out = Vec::new();
    let mut scratch = XScratch::new(false);
    for (f, s_frame) in inp.s_frames_preliminary.iter().enumerate() {
        let chunk = SubjectChunk {
            offset: 0,
            length: s_frame.aa_len,
            overlap: 0,
        };
        let mut init_hsps = Vec::new();
        // The reference finds no hits in a frame shorter than a word.
        if chunk.length >= inp.wordsize as usize {
            for (g, group) in pp.groups.iter().enumerate() {
                if group.lookup.longest_chain > 0 {
                    let unit = pp.x_unit(inp, f, 0, chunk, g);
                    init_hsps.extend(pp.run_unit_batched::<false>(&unit, inp, &mut scratch));
                }
            }
        }
        sort_init_hsps_by_score_ncbi(&mut init_hsps);
        out.push(init_hsps);
    }
    out
}

/// A frame with mask intervals (start, end exclusive, as SEG gives them); the lookup locations
/// are their complement.
fn masked_frame(f: i8, residues: &[u8], masks: &[(usize, usize)]) -> QueryFrame {
    let mut masked = frame(f, residues);
    masked.seg_masks = masks.to_vec();
    masked
}

/// A random frame: now and then empty or shorter than a word, with zero to two mask intervals
/// inside it, now and then masked from end to end.
fn random_masked_frame(rng: &mut Rng, f: i8, max_len: u64, alphabet: &[u8]) -> QueryFrame {
    let len = match rng.below(8) {
        0 => 0,
        1 => 2,
        _ => 3 + rng.below(max_len) as usize,
    };
    let residues = random_residues(rng, len, alphabet);
    let mut masks: Vec<(usize, usize)> = Vec::new();
    if len > 0 {
        match rng.below(5) {
            0 => masks.push((0, len)),
            1 => {}
            _ => {
                let a = rng.below(len as u64) as usize;
                let b = a + 1 + rng.below((len - a) as u64) as usize;
                masks.push((a, b));
                if b + 2 < len && rng.below(2) == 0 {
                    let c = b + 1 + rng.below((len - b - 1) as u64) as usize;
                    let d = c + 1 + rng.below((len - c) as u64) as usize;
                    masks.push((c, d));
                }
            }
        }
    }
    masked_frame(f, &residues, &masks)
}

/// The direct groups of `queries` (`lut_direct_groups`, then the stub) against `build_groups`
/// of the reference builder's table, with `lut_shadow_compare` (every
/// field of every group table, `q_min`, `q_max`, the stub's header fields and `longest_chain`,
/// the contexts). Returns the number of groups and of groups with entries.
fn lut_check(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> (usize, usize) {
    lut_check_for(queries, threshold, pool, ScoringMatrix::Blosum62, 3)
}

/// `lut_check` with the scores of `matrix` and words of `word_length` letters (`-matrix`,
/// `-word_size`): the reference builder's table of `build_ncbi_lookup_for_profile`, the
/// direct groups from `XPairLutQuery::x_prepare_for`.
fn lut_check_for(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    pool: &crate::utils::threading::SearchPool<'_>,
    matrix: ScoringMatrix,
    word_length: usize,
) -> (usize, usize) {
    let params = lookup_protein_params_ungapped(matrix);
    let (full, contexts) = crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile(
        queries,
        threshold,
        &params,
        true,
        matrix,
        word_length,
    );
    assert_eq!(full.word_length, word_length as i32);
    let reference = build_groups(&full, &contexts, X_MAX_GROUPS);
    let q = XPairLutQuery::x_prepare_for(queries, threshold, true, matrix, word_length)
        .expect("direct builder applies");
    assert_eq!(q.word_length(), word_length);
    // On the pool when the caller has one, else in order on this thread; words of 4 letters by
    // the joint builder (`lut_wide_groups`), as `lut_direct` dispatches them.
    let direct = if crate::algorithm::tblastx::lookup::x_pair_lut::x_direct_covers(word_length) {
        lut_direct_groups(&q, X_MAX_GROUPS, pool, true)
    } else {
        lut_wide_groups(&q, &group_starts(q.contexts(), X_MAX_GROUPS), pool)
    }
    .expect("every location inside one group");
    let stub = q.x_into_stub(direct.iter().map(|g| &g.lookup));
    let lut = XLutBatch {
        groups: direct,
        sliced: None,
        stub: Some(stub),
    };
    lut_shadow_compare(&lut, &reference, None, &full, &contexts);
    (
        reference.len(),
        reference
            .iter()
            .filter(|g| g.lookup.longest_chain > 0)
            .count(),
    )
}

const FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456,461-465
#[test]
fn direct_groups_equal_build_groups_for_one_query() {
    let mut rng = Rng(0x0dd5_eed5_1234_5678);
    let alphabets = [aa("ACDEKW"), aa("AW"), aa("ACDEFGHIKLMNPQRSTVWY")];
    let serial = crate::utils::threading::SearchPool::x_serial();
    let mut with_entries = 0usize;
    let mut empty_groups = 0usize;
    for case_no in 0..90 {
        let alphabet = &alphabets[case_no % alphabets.len()];
        // 1 to 6 contexts: one group per context.
        let n_ctx = 1 + case_no % 6;
        let frames: Vec<QueryFrame> = (0..n_ctx)
            .map(|c| random_masked_frame(&mut rng, FRAMES[c], 300, alphabet))
            .collect();
        for threshold in [0, 11, 13] {
            let (groups, filled) = lut_check(std::slice::from_ref(&frames), threshold, &serial);
            assert_eq!(groups, n_ctx);
            with_entries += filled;
            empty_groups += groups - filled;
        }
    }
    assert!(with_entries > 300, "groups with entries {with_entries}");
    assert!(empty_groups > 100, "empty groups {empty_groups}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-590
// A threshold of INT_MIN (the conversion of a threshold beyond Int4): every word is a
// neighbour of every query word, so every cell of a group holds every offset of the group.
#[test]
fn direct_groups_equal_build_groups_at_threshold_int_min() {
    let mut rng = Rng(0x7357_0000_0000_0001);
    let alphabet = aa("ACDEFGHIKLMNPQRSTVWY");
    let serial = crate::utils::threading::SearchPool::x_serial();
    for n_ctx in [1usize, 4, 6] {
        let frames: Vec<QueryFrame> = (0..n_ctx)
            .map(|c| random_masked_frame(&mut rng, FRAMES[c], 9, &alphabet))
            .collect();
        let (groups, _) = lut_check(std::slice::from_ref(&frames), i32::MIN, &serial);
        assert_eq!(groups, n_ctx);
    }
    // Two queries: twelve contexts in six groups.
    let queries: Vec<Vec<QueryFrame>> = (0..2)
        .map(|_| {
            FRAMES
                .iter()
                .map(|&f| random_masked_frame(&mut rng, f, 6, &alphabet))
                .collect()
        })
        .collect();
    let (groups, _) = lut_check(&queries, i32::MIN, &serial);
    assert!((2..=X_MAX_GROUPS).contains(&groups), "groups {groups}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1047
// Several queries make more than X_MAX_GROUPS contexts: a group then holds several contexts,
// and its direct table is built from the lookup locations of all of them.
#[test]
fn direct_groups_equal_build_groups_for_several_queries() {
    let mut rng = Rng(0x0bad_cafe_dead_beef);
    let alphabets = [aa("ACDEKW"), aa("ACDEFGHIKLMNPQRSTVWY")];
    let serial = crate::utils::threading::SearchPool::x_serial();
    let mut six_groups = 0usize;
    for case_no in 0..24 {
        let alphabet = &alphabets[case_no % alphabets.len()];
        let n_queries = 2 + case_no % 5;
        let queries: Vec<Vec<QueryFrame>> = (0..n_queries)
            .map(|_| {
                FRAMES
                    .iter()
                    .map(|&f| random_masked_frame(&mut rng, f, 120, alphabet))
                    .collect()
            })
            .collect();
        for threshold in [0, 11, 13] {
            // 12 to 36 contexts in at most X_MAX_GROUPS groups (`group_starts` balances them
            // by length and can make fewer), so some group holds several contexts.
            let (groups, _) = lut_check(&queries, threshold, &serial);
            assert!((1..=X_MAX_GROUPS).contains(&groups), "groups {groups}");
            six_groups += usize::from(groups == X_MAX_GROUPS);
        }
    }
    assert!(six_groups > 0, "cases with six groups {six_groups}");
    // Four queries with frames of one length: six groups of four contexts.
    let even: Vec<Vec<QueryFrame>> = (0..4)
        .map(|_| {
            FRAMES
                .iter()
                .map(|&f| masked_frame(f, &random_residues(&mut rng, 60, &alphabets[1]), &[]))
                .collect()
        })
        .collect();
    assert_eq!(lut_check(&even, 13, &serial), (X_MAX_GROUPS, X_MAX_GROUPS));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:96-100
// A context masked from end to end has no lookup location, an empty frame has no letter and a
// two-letter frame no word: their groups have no entry, `q_min` = `q_max` = 0, and the units
// skip them (`longest_chain` 0).
#[test]
fn direct_groups_of_contexts_without_words_are_empty() {
    let mut rng = Rng(0x00c0_ffee_0000_0042);
    let alphabet = aa("ACDEKLMNPQRSTVWY");
    let serial = crate::utils::threading::SearchPool::x_serial();
    let frames = vec![
        masked_frame(1, &random_residues(&mut rng, 80, &alphabet), &[(10, 20)]),
        masked_frame(2, &random_residues(&mut rng, 70, &alphabet), &[(0, 70)]),
        masked_frame(3, &[], &[]),
        masked_frame(-1, &random_residues(&mut rng, 2, &alphabet), &[]),
        masked_frame(-2, &random_residues(&mut rng, 60, &alphabet), &[]),
        masked_frame(
            -3,
            &random_residues(&mut rng, 50, &alphabet),
            &[(0, 3), (47, 50)],
        ),
    ];
    for threshold in [0, 11, 13] {
        assert_eq!(
            lut_check(std::slice::from_ref(&frames), threshold, &serial),
            (6, 3)
        );
    }
    let q = XPairLutQuery::x_prepare(std::slice::from_ref(&frames), 13, true)
        .expect("direct builder applies");
    let groups = lut_direct_groups(&q, X_MAX_GROUPS, &serial, false)
        .expect("every location inside one group");
    for g in [1usize, 2, 3] {
        assert_eq!(groups[g].lookup.longest_chain, 0);
        assert_eq!((groups[g].q_min, groups[g].q_max), (0, 0));
        assert!(groups[g].lookup.overflow.is_empty());
    }
    // The stub has the full table's longest chain although it has no entry.
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (full, _) = build_ncbi_lookup(std::slice::from_ref(&frames), 13, &params, true);
    let (stub, contexts) = q.x_into_stub(groups.iter().map(|g| &g.lookup));
    assert!(full.longest_chain > 0);
    assert_eq!(stub.longest_chain, full.longest_chain);
    assert_eq!(contexts.len(), 6);
    assert!(stub.pv.iter().all(|&word| word == 0) && stub.overflow.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
// The groups built on a search pool (each build on a memo of neighbour lists of
// `X_TBX_MEMOS`) are the groups built in order on one thread.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
#[test]
fn direct_groups_on_a_search_pool_equal_build_groups() {
    let mut rng = Rng(0x9001_9002_9003_9004);
    let alphabet = aa("ACDEFGHIKLMNPQRSTVWY");
    let one: Vec<Vec<QueryFrame>> = vec![FRAMES
        .iter()
        .map(|&f| random_masked_frame(&mut rng, f, 2000, &alphabet))
        .collect()];
    let many: Vec<Vec<QueryFrame>> = (0..4)
        .map(|_| {
            FRAMES
                .iter()
                .map(|&f| random_masked_frame(&mut rng, f, 400, &alphabet))
                .collect()
        })
        .collect();
    crate::utils::threading::with_search_pool(4, "tblastx", |pool| {
        assert!(pool.enabled());
        for threshold in [0, 11, 13] {
            // The second round reads the lists memoised in the first.
            for _round in 0..2 {
                assert_eq!(lut_check(&one, threshold, pool).0, 6);
                let (groups, _) = lut_check(&many, threshold, pool);
                assert!((2..=X_MAX_GROUPS).contains(&groups), "groups {groups}");
            }
        }
        Ok(())
    })
    .expect("search pool");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
// Builds of two keys side by side (two searches in one process) each take a memo of
// `X_TBX_MEMOS` and give it back; a memo kept for the other key is replaced, never read.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
#[test]
fn direct_groups_of_two_keys_at_once_equal_build_groups() {
    let mut rng = Rng(0x9005_9006_9007_9008);
    let alphabet = aa("ACDEFGHIKLMNPQRSTVWY");
    let many: Vec<Vec<QueryFrame>> = (0..4)
        .map(|_| {
            FRAMES
                .iter()
                .map(|&f| random_masked_frame(&mut rng, f, 400, &alphabet))
                .collect()
        })
        .collect();
    crate::utils::threading::with_search_pool(4, "tblastx", |pool| {
        for _round in 0..3 {
            let (a, b) = pool.install(|| {
                rayon::join(|| lut_check(&many, 11, pool), || lut_check(&many, 13, pool))
            });
            assert!((2..=X_MAX_GROUPS).contains(&a.0), "groups {}", a.0);
            assert!((2..=X_MAX_GROUPS).contains(&b.0), "groups {}", b.0);
        }
        Ok(())
    })
    .expect("search pool");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-57
// `decide` is the decision of `prepare`: the gates, then the guard.
#[test]
fn decide_checks_the_gates_and_the_guard() {
    let mut rng = Rng(5);
    let alphabet = aa("ACDEKW");
    let window = 40;
    // One context of 216 letters: the table has 256 cells and the slack is 40.
    let (_, contexts) = lookup_for(&[frame(1, &random_residues(&mut rng, 216, &alphabet))]);
    let query_length = query_length_of(&contexts);
    let diag_array_size = diag_array_size_for(query_length, window);
    assert_eq!((query_length, diag_array_size), (216, 256));
    let short = [FastaRecord::new("", b"s", &[b'A'; 120])];
    let long = [
        FastaRecord::new("", b"s", &[b'A'; 30]),
        FastaRecord::new("", b"l", &[b'A'; 123]),
    ];
    let none: [FastaRecord; 0] = [];
    let decision_band =
        |subjects: &[FastaRecord], reference_loop: bool, chunk_switches: bool, band: bool| {
            decide(&XGates {
                contexts: &contexts,
                query_length,
                diag_array_size,
                subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(subjects),
                reference_loop,
                chunk_switches,
                threads: 1,
                cost_gate: false,
                n1_gate: false,
                threshold: 13,
                band,
                unit_mask: false,
            })
        };
    let decision = |subjects: &[FastaRecord], reference_loop: bool, chunk_switches: bool| {
        decision_band(subjects, reference_loop, chunk_switches, false)
    };
    assert_eq!(decision(&short, false, false), XDecision::Units);
    assert_eq!(decision(&none, false, false), XDecision::Units);
    assert_eq!(decision(&long, false, false), XDecision::Guard);
    assert_eq!(decision(&short, true, false), XDecision::Gate);
    assert_eq!(decision(&short, false, true), XDecision::Gate);
    // The gates come before the guard, as in the counters of LOSAT_X_STATS.
    assert_eq!(decision(&long, true, false), XDecision::Gate);
    // LOSAT_X_TBXBAND: the first guard is not applied; the gates still are.
    assert_eq!(decision_band(&long, false, false, true), XDecision::Units);
    assert_eq!(decision_band(&short, false, false, true), XDecision::Units);
    assert_eq!(decision_band(&long, true, false, true), XDecision::Gate);
    assert_eq!(decision_band(&short, false, true, true), XDecision::Gate);
    // The memory rule: a subject chunk longer than the table's 256 cells keeps the guard.
    let longer = [FastaRecord::new("", b"l", &[b'A'; 3 * 257])];
    let fits = [FastaRecord::new("", b"l", &[b'A'; 3 * 256])];
    assert_eq!(decision_band(&longer, false, false, true), XDecision::Guard);
    assert_eq!(decision_band(&fits, false, false, true), XDecision::Units);
    assert!(x_band_fits(256, 256) && !x_band_fits(257, 256));
    // LOSAT_X_TBXUNITMASK: the memory rule is not applied where LOSAT_X_TBXQSLICE cuts no group
    // (small groups at one thread); the gates and the second guard still are.
    let decision_mask = |subjects: &[FastaRecord], reference_loop: bool| {
        decide(&XGates {
            contexts: &contexts,
            query_length,
            diag_array_size,
            subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(subjects),
            reference_loop,
            chunk_switches: false,
            threads: 1,
            cost_gate: false,
            n1_gate: false,
            threshold: 13,
            band: true,
            unit_mask: true,
        })
    };
    assert_eq!(decision_mask(&longer, false), XDecision::Units);
    assert_eq!(decision_mask(&fits, false), XDecision::Units);
    assert_eq!(decision_mask(&longer, true), XDecision::Gate);
    // LOSAT_X_PAIRPARGATE still decides after the guard: this table of 256 cells at one thread.
    let gated = |band: bool| {
        decide(&XGates {
            contexts: &contexts,
            query_length,
            diag_array_size,
            subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(&long),
            reference_loop: false,
            chunk_switches: false,
            threads: 1,
            cost_gate: true,
            n1_gate: false,
            threshold: 13,
            band,
            unit_mask: false,
        })
    };
    assert_eq!(gated(true), XDecision::Cost);
    assert_eq!(gated(false), XDecision::Guard);
    // The second guard (`last_hit_fits`): with a window just under 2^29 the slack is almost
    // 2^29, so the first guard holds for the short subject and the units still do not run.
    let wide = diag_array_size_for(query_length, 536_870_000);
    assert_eq!(wide, 1 << 29);
    assert!(guard(216, 40, query_length, wide));
    let gates = XGates {
        contexts: &contexts,
        query_length,
        diag_array_size: wide,
        subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(&short),
        reference_loop: false,
        chunk_switches: false,
        threads: 1,
        cost_gate: false,
        n1_gate: false,
        threshold: 13,
        band: false,
        unit_mask: false,
    };
    assert_eq!(decide(&gates), XDecision::Guard);
    // LOSAT_X_TBXBAND does not lift the second guard.
    assert_eq!(
        decide(&XGates {
            band: true,
            ..gates
        }),
        XDecision::Guard
    );
}

// No NCBI counterpart: LOSAT_X_PAIRPARGATE chooses between two LOSAT paths with the same result.
// The cost gate comes after the gates and the guard, keeps the units below
// the default threshold, for a reference table of 2^18 cells or more, and at several threads
// for a subject that holds more than one and a half thread's share of the bases.
#[test]
fn the_cost_gate_keeps_the_units_where_they_gain() {
    let mut rng = Rng(7);
    let alphabet = aa("ACDEKW");
    let window = 40;
    let (_, contexts) = lookup_for(&[frame(1, &random_residues(&mut rng, 216, &alphabet))]);
    let query_length = query_length_of(&contexts);
    let diag_array_size = diag_array_size_for(query_length, window);
    assert_eq!(diag_array_size, 256);
    let one = [FastaRecord::new("", b"s", &[b'A'; 120])];
    // 120 of 200 bases in one subject: 2 * 120 * t > 3 * 200 from t = 3.
    let mixed = [
        FastaRecord::new("", b"s", &[b'A'; 120]),
        FastaRecord::new("", b"t", &[b'A'; 40]),
        FastaRecord::new("", b"u", &[b'A'; 40]),
    ];
    let none: [FastaRecord; 0] = [];
    let long = [
        FastaRecord::new("", b"s", &[b'A'; 30]),
        FastaRecord::new("", b"l", &[b'A'; 123]),
    ];
    let decision_at = |threshold: i32,
                       diag_array_size: i32,
                       subjects: &[FastaRecord],
                       reference_loop: bool,
                       threads: usize,
                       cost_gate: bool| {
        decide(&XGates {
            contexts: &contexts,
            query_length,
            diag_array_size,
            subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(subjects),
            reference_loop,
            chunk_switches: false,
            threads,
            cost_gate,
            n1_gate: false,
            threshold,
            band: false,
            unit_mask: false,
        })
    };
    let decision = |diag_array_size: i32,
                    subjects: &[FastaRecord],
                    reference_loop: bool,
                    threads: usize,
                    cost_gate: bool| {
        decision_at(
            13,
            diag_array_size,
            subjects,
            reference_loop,
            threads,
            cost_gate,
        )
    };
    // Without the gate nothing changes.
    for t in [1, 2, 8] {
        assert_eq!(decision(256, &one, false, t, false), XDecision::Units);
    }
    // One thread: the table size alone.
    assert_eq!(decision(256, &one, false, 1, true), XDecision::Cost);
    assert_eq!(decision(1 << 17, &one, false, 1, true), XDecision::Cost);
    assert_eq!(decision(1 << 18, &one, false, 1, true), XDecision::Units);
    assert_eq!(decision(256, &none, false, 1, true), XDecision::Cost);
    // Several threads: one subject is all of the bases.
    assert_eq!(decision(256, &one, false, 2, true), XDecision::Units);
    assert_eq!(decision(256, &none, false, 2, true), XDecision::Cost);
    assert_eq!(decision(256, &mixed, false, 2, true), XDecision::Cost);
    assert_eq!(decision(256, &mixed, false, 3, true), XDecision::Units);
    assert_eq!(decision(256, &mixed, false, 8, true), XDecision::Units);
    assert_eq!(decision(1 << 18, &mixed, false, 2, true), XDecision::Units);
    // Below the default threshold (and INT_MIN, `-threshold +inf`) the units stay.
    assert_eq!(decision_at(12, 256, &one, false, 1, true), XDecision::Units);
    assert_eq!(
        decision_at(i32::MIN, 256, &one, false, 1, true),
        XDecision::Units
    );
    assert_eq!(decision_at(14, 256, &one, false, 1, true), XDecision::Cost);
    // The gates and the guard come first.
    assert_eq!(decision(256, &one, true, 1, true), XDecision::Gate);
    assert_eq!(decision(256, &long, false, 1, true), XDecision::Guard);
}

// No NCBI counterpart: LOSAT_X_TBXFEWGROUPS's groups: two units per thread on the longest
// subject's share of the bases, at most six groups, and only up to the letters' line.
#[test]
fn a_batch_takes_the_groups_the_pool_needs() {
    let subjects = |lengths: &[usize]| -> Vec<FastaRecord> {
        lengths
            .iter()
            .map(|&length| FastaRecord::new("", b"s", &vec![b'A'; length]))
            .collect()
    };
    let wanted = |threads: usize, records: &[FastaRecord]| {
        x_groups_wanted(
            threads,
            crate::algorithm::tblastx::subjects::Subjects::Fasta(records),
        )
    };
    let one = subjects(&[900]);
    assert_eq!(wanted(1, &one), 1);
    assert_eq!(wanted(2, &one), 1);
    assert_eq!(wanted(3, &one), 1);
    assert_eq!(wanted(4, &one), 2);
    assert_eq!(wanted(8, &one), 3);
    assert_eq!(wanted(16, &one), 6);
    assert_eq!(wanted(48, &one), 6);
    let equal = subjects(&[900; 40]);
    assert_eq!(wanted(8, &equal), 1);
    assert_eq!(wanted(48, &equal), 1);
    // A chromosome among plasmids needs the units of a subject alone.
    let chromosome = subjects(&[50_000, 500, 500, 500, 500, 500]);
    assert_eq!(wanted(16, &chromosome), 6);
    assert_eq!(wanted(8, &chromosome), 3);
    let mut rng = Rng(5);
    let alphabet = aa("ACDEKW");
    let frames: Vec<QueryFrame> = (0..6)
        .map(|k| frame(1 + k, &random_residues(&mut rng, 30, &alphabet)))
        .collect();
    let (_, contexts) = lookup_for(&frames);
    assert_eq!(group_starts(&contexts, 1), vec![0]);
    assert_eq!(group_starts(&contexts, 3), vec![0, 2, 4]);
    assert_eq!(
        group_starts(&contexts, X_MAX_GROUPS),
        vec![0, 1, 2, 3, 4, 5]
    );
    assert_eq!(x_takes_few_groups(&contexts, 3), Some(true));
    assert_eq!(x_takes_few_groups(&contexts, X_MAX_GROUPS), None);
}

// No NCBI counterpart: LOSAT_X_PAIRPARN1 moves the one-thread line of LOSAT_X_PAIRPARGATE from
// 2^18 cells to batches of 5,000 letters, and LOSAT_X_TBXFEWGROUPS takes the same line at several
// threads; the threshold rule and the gates before them are unchanged.
#[test]
fn the_one_thread_line_keeps_the_units_from_5000_letters() {
    let mut rng = Rng(11);
    let alphabet = aa("ACDEKW");
    let (_, contexts) = lookup_for(&[frame(1, &random_residues(&mut rng, 216, &alphabet))]);
    let one = [FastaRecord::new("", b"s", &[b'A'; 120])];
    // 10 subjects of 40 bases: no subject holds one and a half thread's share at 2 or 8 threads.
    let spread: Vec<FastaRecord> = (0..10)
        .map(|_| FastaRecord::new("", b"s", &[b'A'; 40]))
        .collect();
    // The batch's length in letters (the table follows it, as at the build site).
    let decision = |query_length: i32,
                    subjects: &[FastaRecord],
                    threads: usize,
                    threshold: i32,
                    n1_gate: bool,
                    reference_loop: bool| {
        decide(&XGates {
            contexts: &contexts,
            query_length,
            diag_array_size: diag_array_size_for(query_length, 40),
            subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(subjects),
            reference_loop,
            chunk_switches: false,
            threads,
            cost_gate: true,
            n1_gate,
            threshold,
            band: false,
            unit_mask: false,
        })
    };
    assert_eq!(diag_array_size_for(140_000, 40), 1 << 18);
    // One thread: the units from 5,000 letters with the switch, from 2^18 cells without it.
    assert_eq!(decision(4_999, &one, 1, 13, true, false), XDecision::Cost);
    assert_eq!(decision(5_000, &one, 1, 13, true, false), XDecision::Units);
    assert_eq!(
        decision(100_000, &spread, 1, 13, true, false),
        XDecision::Units
    );
    assert_eq!(decision(5_000, &one, 1, 13, false, false), XDecision::Cost);
    assert_eq!(
        decision(100_000, &one, 1, 13, false, false),
        XDecision::Cost
    );
    assert_eq!(
        decision(140_000, &one, 1, 13, false, false),
        XDecision::Units
    );
    // Several threads, LOSAT_X_TBXFEWGROUPS: the same line for spread subjects.
    for threads in [2, 8] {
        assert_eq!(
            decision(4_999, &spread, threads, 13, true, false),
            XDecision::Cost
        );
        assert_eq!(
            decision(5_000, &spread, threads, 13, true, false),
            XDecision::Units
        );
        assert_eq!(
            decision(5_000, &one, threads, 13, true, false),
            XDecision::Units
        );
    }
    // The threshold rule and the gates come first.
    assert_eq!(decision(4_999, &one, 1, 12, true, false), XDecision::Units);
    assert_eq!(decision(5_000, &one, 1, 13, true, true), XDecision::Gate);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:264-266
// The build site gets the stub once.
#[test]
fn the_stub_is_taken_once() {
    let mut rng = Rng(9);
    let alphabet = aa("ACDEKW");
    let serial = crate::utils::threading::SearchPool::x_serial();
    let queries = vec![vec![
        frame(1, &random_residues(&mut rng, 50, &alphabet)),
        frame(2, &random_residues(&mut rng, 40, &alphabet)),
    ]];
    let batch = || {
        let q = XPairLutQuery::x_prepare(&queries, 11, true).expect("direct builder applies");
        let groups = lut_direct_groups(&q, X_MAX_GROUPS, &serial, false)
            .expect("every location inside one group");
        let stub = q.x_into_stub(groups.iter().map(|g| &g.lookup));
        XLutBatch {
            groups,
            sliced: None,
            stub: Some(stub),
        }
    };
    let mut direct = batch();
    let (stub, contexts) = direct.take_stub().expect("the stub of the batch");
    assert_eq!(contexts.len(), 2);
    assert_eq!(stub.num_contexts, 2);
    assert!(direct.take_stub().is_none());
    assert!(direct.stub.is_none());
}

// No NCBI counterpart: grouping of contexts (scheduling only).
#[test]
fn many_contexts_make_six_contiguous_groups() {
    let mut rng = Rng(3);
    let alphabet = aa("ACDEKW");
    let queries: Vec<Vec<QueryFrame>> = (0..4)
        .map(|_| {
            [1i8, 2, 3, -1, -2, -3]
                .iter()
                .map(|&f| frame(f, &random_residues(&mut rng, 20, &alphabet)))
                .collect()
        })
        .collect();
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (_lookup, contexts) = build_ncbi_lookup(&queries, 11, &params, true);
    let starts = group_starts(&contexts, X_MAX_GROUPS);
    assert_eq!(starts.len(), X_MAX_GROUPS);
    assert_eq!(starts[0], 0);
    assert!(starts.windows(2).all(|w| w[0] < w[1]));
    assert!(*starts.last().unwrap() < contexts.len());
}

/// One random case of the LOSAT_X_SEEDCOMPACT tests: one query of one to six contexts (now
/// and then two or three queries, so that a group holds several contexts), its groups, and one
/// to six subject frames of `s_min` to `s_min + s_spread - 1` letters (now and then shorter
/// than a word).
struct CompactCase {
    pp: XPairPar,
    contexts: Vec<QueryContext>,
    s_frames: Vec<QueryFrame>,
    cutoff_scores: Vec<i32>,
    x_dropoff_per_context: Vec<i32>,
    failure: crate::algorithm::tblastx::subjects::XReadFailure,
}

impl CompactCase {
    fn random(rng: &mut Rng, alphabet: &[u8], s_min: usize, s_spread: u64) -> Self {
        let n_queries = if rng.below(4) == 0 {
            2 + rng.below(2) as usize
        } else {
            1
        };
        let mut queries: Vec<Vec<QueryFrame>> = Vec::new();
        for _ in 0..n_queries {
            let n_ctx = 1 + rng.below(6) as usize;
            let mut frames = Vec::new();
            for &f in &FRAMES[..n_ctx] {
                let len = 3 + rng.below(60) as usize;
                frames.push(frame(f, &random_residues(rng, len, alphabet)));
            }
            queries.push(frames);
        }
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (lookup, contexts) = build_ncbi_lookup(&queries, 11, &params, true);
        let n_s = 1 + rng.below(6) as usize;
        let mut s_frames = Vec::new();
        for &f in &FRAMES[..n_s] {
            let len = if rng.below(10) == 0 {
                rng.below(3) as usize
            } else {
                s_min + rng.below(s_spread) as usize
            };
            s_frames.push(frame(f, &random_residues(rng, len, alphabet)));
        }
        Self {
            pp: XPairPar {
                groups: build_groups(&lookup, &contexts, X_MAX_GROUPS),
                stub: false,
                sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
                sc: Blosum62Score,
                sliced: None,
                mask_rule: None,
                extents: Vec::new(),
                ctx_of: x_ctx_table(&lookup, true),
            },
            cutoff_scores: vec![12; contexts.len()],
            x_dropoff_per_context: vec![16; contexts.len()],
            contexts,
            s_frames,
            failure: crate::algorithm::tblastx::subjects::XReadFailure::new(),
        }
    }

    fn inp(&self) -> XSubject<'_> {
        XSubject {
            s_frames_preliminary: &self.s_frames,
            s_frames: &self.s_frames,
            contexts: &self.contexts,
            cutoff_scores: &self.cutoff_scores,
            x_dropoff_per_context: &self.x_dropoff_per_context,
            window: 40,
            wordsize: 3,
            s_idx: 7,
            s_len: 4_321,
            failure: &self.failure,
            passed: None,
            subject_besthit: false,
        }
    }
}

/// The chunks of every subject frame and the units of every (frame, chunk), listed with the
/// statements of `run_subject_units` for chunks of at most `max_dbseq_len` letters
/// (`seed_compact_test_listing_is_run_subject_units` checks it at MAX_DBSEQ_LEN).
fn list_units(
    pp: &XPairPar,
    inp: &XSubject<'_>,
    max_dbseq_len: usize,
) -> (Vec<Vec<SubjectChunk>>, Vec<XUnit>) {
    let mut frame_chunks = Vec::new();
    let mut units = Vec::new();
    for (f, s_frame) in inp.s_frames_preliminary.iter().enumerate() {
        let mut split_state = SubjectSplitState::new(s_frame.aa_len);
        let mut chunks = Vec::new();
        while let SubjectChunkStatus::Ok(chunk) =
            split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP)
        {
            chunks.push(chunk);
        }
        for (k, chunk) in chunks.iter().enumerate() {
            if chunk.length < inp.wordsize as usize {
                continue;
            }
            let first = pp.x_unit(inp, f, k, *chunk, 0);
            for (g, group) in pp.groups.iter().enumerate() {
                if group.lookup.longest_chain > 0 {
                    units.push(XUnit { g, ..first });
                }
            }
        }
        frame_chunks.push(chunks);
    }
    (frame_chunks, units)
}

/// What `compact_check` saw: units, units without a saved init HSP, compact records, frames
/// whose list got hits from two or more chunks (the chunk merge ran), and hits of the subject.
#[derive(Default)]
struct CompactSeen {
    units: usize,
    empty_units: usize,
    records: usize,
    merged_frames: usize,
    hits: usize,
}

impl CompactSeen {
    fn add(&mut self, other: CompactSeen) {
        self.units += other.units;
        self.empty_units += other.empty_units;
        self.records += other.records;
        self.merged_frames += other.merged_frames;
        self.hits += other.hits;
    }
}

/// The units of `case` (chunks of at most `max_dbseq_len` letters) run on `pool` both ways:
/// every unit's compact records have no spare room and expand to today's list of that unit,
/// element for element, and `merge_subject_listed` of the compact records gives today's
/// `combined_ungapped_hits`, every field (`shadow_compare`), as does
/// the merge of LOSAT_X_TBXKMERGE of either form, sorted as the units' tasks sort them.
fn compact_check(
    case: &CompactCase,
    pool: &crate::utils::threading::SearchPool<'_>,
    max_dbseq_len: usize,
) -> CompactSeen {
    let inp = case.inp();
    let (frame_chunks, units) = list_units(&case.pp, &inp, max_dbseq_len);
    let full = case
        .pp
        .run_units(&units, &inp, pool, XPairPar::run_unit::<false>, false);
    let compact = case.pp.run_units(
        &units,
        &inp,
        pool,
        XPairPar::run_unit_compact::<false>,
        false,
    );
    assert_eq!((full.len(), compact.len()), (units.len(), units.len()));
    let mut seen = CompactSeen {
        units: units.len(),
        ..CompactSeen::default()
    };
    for (unit, (list, records)) in units.iter().zip(full.iter().zip(&compact)) {
        let at = format!("unit f {} k {} g {}", unit.f, unit.k, unit.g);
        // No spare room, except in a copy that LOSAT_X_SEEDUNMAP gave the mapped capacity.
        let mapped_copy = records.len() * std::mem::size_of::<XUnitHsp>() >= X_MAPPED_MIN
            && records.capacity() == x_mapped_len::<XUnitHsp>(records.len());
        assert!(
            records.capacity() == records.len() || mapped_copy,
            "{at}: spare room"
        );
        let mut expanded = Vec::new();
        x_expand_unit(records, unit.f, &inp, &mut expanded);
        assert_eq!(&expanded, list, "{at}");
        seen.empty_units += usize::from(list.is_empty());
        seen.records += records.len();
    }
    for (f, chunks) in frame_chunks.iter().enumerate() {
        let chunks_with_hits = (0..chunks.len())
            .filter(|&k| {
                units
                    .iter()
                    .zip(&full)
                    .any(|(u, list)| u.f == f && u.k == k && !list.is_empty())
            })
            .count();
        seen.merged_frames += usize::from(chunks_with_hits >= 2);
    }
    // LOSAT_X_TBXKMERGE (x_kmerge.rs): the units' lists sorted as their tasks sort them, merged
    // by the switch's merge; both forms give today's list.
    let mut scratch = XScratch::new(false);
    let mut sorted_full = full.clone();
    sorted_full
        .iter_mut()
        .for_each(|list| list.x_sort_in(&mut scratch));
    let mut sorted_compact = compact.clone();
    sorted_compact
        .iter_mut()
        .for_each(|records| records.x_sort_in(&mut scratch));
    let today = XPairPar::<Blosum62Score>::merge_subject_listed(
        XSubjectUnits {
            frame_chunks: frame_chunks.clone(),
            units: units.clone(),
            outputs: XUnitOutputs::Full(full),
        },
        &inp,
    );
    let with_switch = XPairPar::<Blosum62Score>::merge_subject_listed(
        XSubjectUnits {
            frame_chunks: frame_chunks.clone(),
            units: units.clone(),
            outputs: XUnitOutputs::Compact(compact),
        },
        &inp,
    );
    shadow_compare(inp.s_idx, &with_switch, &today);
    for outputs in [
        XUnitOutputs::Full(sorted_full),
        XUnitOutputs::Compact(sorted_compact),
    ] {
        let merged = x_kmerge::merge_subject(
            XSubjectUnits {
                frame_chunks: frame_chunks.clone(),
                units: units.clone(),
                outputs,
            },
            &inp,
        );
        shadow_compare(inp.s_idx, &merged, &today);
    }
    seen.hits = today.len();
    seen
}

const COMPACT_ALPHABETS: [&str; 3] = ["ACDEKW", "AW", "ACDEFGHIKLMNPQRSTVWY"];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:367-372
// LOSAT_X_SEEDCOMPACT on random units, one chunk per frame, serially (one scratch for every
// unit of a case): each unit's records are its init HSPs of today, and the merge is today's.
#[test]
fn seed_compact_unit_records_expand_to_todays_lists() {
    let mut rng = Rng(0x5eed_c0de_0000_0001);
    let serial = crate::utils::threading::SearchPool::x_serial();
    let mut seen = CompactSeen::default();
    for case_no in 0..150 {
        let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
        let case = CompactCase::random(&mut rng, &alphabet, 3, 400);
        seen.add(compact_check(&case, &serial, MAX_DBSEQ_LEN));
    }
    assert!(seen.units > 1_000, "units {}", seen.units);
    assert!(seen.empty_units > 50, "empty units {}", seen.empty_units);
    assert!(seen.records > 2_000, "records {}", seen.records);
    assert!(seen.hits > 1_000, "hits {}", seen.hits);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2876-2881
// Subjects of several chunks per frame (chunks of 103, 160 and 300 letters with the overlap
// of 100: the first chunk's list moves into the empty frame list, the later ones go through
// `merge_tblastx_subject_chunk_hits`): the compact merge gives today's list, every field.
#[test]
fn seed_compact_merge_with_several_chunks_per_frame_is_todays_merge() {
    let mut rng = Rng(0x5eed_c0de_0000_0002);
    let serial = crate::utils::threading::SearchPool::x_serial();
    let mut seen = CompactSeen::default();
    for case_no in 0..36 {
        let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
        let case = CompactCase::random(&mut rng, &alphabet, 200, 500);
        let max_dbseq_len = [DBSEQ_CHUNK_OVERLAP + 3, 160, 300][case_no % 3];
        seen.add(compact_check(&case, &serial, max_dbseq_len));
    }
    assert!(
        seen.merged_frames > 20,
        "merged frames {}",
        seen.merged_frames
    );
    assert!(seen.hits > 1_000, "hits {}", seen.hits);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1043
// On search pools of 1, 2 and 8 threads (one scratch vector per rayon job, kept across its
// units): the same records, lists and merge as today, with one and with several chunks per
// frame.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
#[test]
fn seed_compact_lists_on_search_pools_are_todays() {
    for threads in [1usize, 2, 8] {
        crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
            assert_eq!(pool.enabled(), threads > 1);
            let mut rng = Rng(0x5eed_c0de_0000_0003);
            let mut seen = CompactSeen::default();
            for case_no in 0..30 {
                let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
                let (case, max_dbseq_len) = if case_no % 2 == 0 {
                    (
                        CompactCase::random(&mut rng, &alphabet, 3, 400),
                        MAX_DBSEQ_LEN,
                    )
                } else {
                    (CompactCase::random(&mut rng, &alphabet, 200, 400), 160)
                };
                seen.add(compact_check(&case, pool, max_dbseq_len));
            }
            assert!(
                seen.records > 500,
                "threads {threads}: records {}",
                seen.records
            );
            assert!(
                seen.merged_frames > 5,
                "threads {threads}: merged frames {}",
                seen.merged_frames
            );
            Ok(())
        })
        .expect("search pool");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2371-2375
// One record, through the sink and back, for context indices below and above 16 bits (a
// batch of 70,001 contexts with distinct query indices, frames, lengths and offsets); no
// record expands to nothing. The record is 28 bytes.
#[test]
fn seed_compact_one_record_and_context_indices_above_16_bits_round_trip() {
    assert_eq!(std::mem::size_of::<XUnitHsp>(), 28);
    let (_, template) = lookup_for(&[frame(1, &aa("ACDEKWACDEKW"))]);
    let contexts: Vec<QueryContext> = (0..70_001usize)
        .map(|c| QueryContext {
            q_idx: (c / 6) as u32,
            frame: FRAMES[c % 6],
            orig_len: 100 + c,
            frame_base: (c * 13) as i32,
            ..template[0].clone()
        })
        .collect();
    let s_frames = vec![
        frame(2, &aa("ACDEKWACDEKW")),
        frame(-3, &aa("WKEDCAWKEDCA")),
    ];
    let inp = XSubject {
        s_frames_preliminary: &s_frames,
        s_frames: &s_frames,
        contexts: &contexts,
        cutoff_scores: &[],
        x_dropoff_per_context: &[],
        window: 40,
        wordsize: 3,
        s_idx: 3,
        s_len: 38,
        failure: &crate::algorithm::tblastx::subjects::XReadFailure::new(),
        passed: None,
        subject_besthit: false,
    };
    for (ctx_idx, f) in [(0usize, 0usize), (65_535, 1), (65_536, 0), (70_000, 1)] {
        let ctx = &contexts[ctx_idx];
        // The InitHSP `x_hit` pushes for hsp_q 2, hsp_s 1, hsp_len 5 (its statements).
        let (hsp_q, hsp_s, hsp_len) = (2, 1, 5);
        let hsp = InitHSP {
            q_start_absolute: ctx.frame_base + hsp_q,
            q_end_absolute: ctx.frame_base + (hsp_q + hsp_len),
            s_start: hsp_s,
            s_end: hsp_s + hsp_len,
            q_seed_absolute: ctx.frame_base + 4,
            s_seed: 3,
            score: 31,
            ctx_idx,
            s_f_idx: f,
            q_idx: ctx.q_idx,
            s_idx: inp.s_idx as u32,
            q_frame: ctx.frame,
            s_frame: s_frames[f].frame,
            q_orig_len: ctx.orig_len,
            s_orig_len: inp.s_len,
        };
        let mut records: Vec<XUnitHsp> = Vec::new();
        XInitSink::push(&mut records, hsp);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].ctx_idx as usize, ctx_idx);
        let mut expanded = Vec::new();
        x_expand_unit(&records, f, &inp, &mut expanded);
        assert_eq!(expanded, vec![hsp], "context {ctx_idx}");
        // Today's sink keeps the InitHSP itself.
        let mut today: Vec<InitHSP> = Vec::new();
        XInitSink::push(&mut today, hsp);
        assert_eq!(today, vec![hsp]);
    }
    let mut untouched: Vec<InitHSP> = Vec::new();
    x_expand_unit(&[], 1, &inp, &mut untouched);
    assert!(untouched.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:479-491
// The listing of the tests (`list_units`) is the listing of `run_subject_units` at
// MAX_DBSEQ_LEN, and `run_subject_units` gives the form the switch of this process selects
// (today's lists without LOSAT_X_SEEDCOMPACT), with the lists of `run_units` (sorted by the
// units where LOSAT_X_TBXKMERGE is on, x_kmerge.rs item 1).
#[test]
fn seed_compact_test_listing_is_run_subject_units() {
    let mut rng = Rng(0x5eed_c0de_0000_0004);
    let serial = crate::utils::threading::SearchPool::x_serial();
    for case_no in 0..20 {
        let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
        let case = CompactCase::random(&mut rng, &alphabet, 3, 300);
        let inp = case.inp();
        let x_units = case.pp.run_subject_units(&inp, &serial);
        let (frame_chunks, units) = list_units(&case.pp, &inp, tblastx_max_dbseq_len_for_run());
        assert_eq!(x_units.frame_chunks, frame_chunks);
        let key = |u: &XUnit| (u.f, u.k, u.chunk, u.g);
        assert_eq!(
            x_units.units.iter().map(key).collect::<Vec<_>>(),
            units.iter().map(key).collect::<Vec<_>>()
        );
        match (seed_compact(), x_units.outputs) {
            (false, XUnitOutputs::Full(outputs)) => {
                assert_eq!(
                    outputs,
                    case.pp.run_units(
                        &units,
                        &inp,
                        &serial,
                        XPairPar::run_unit::<false>,
                        x_kmerge::on()
                    )
                );
            }
            (true, XUnitOutputs::Compact(outputs)) => {
                assert_eq!(
                    outputs,
                    case.pp.run_units(
                        &units,
                        &inp,
                        &serial,
                        XPairPar::run_unit_compact::<false>,
                        x_kmerge::on()
                    )
                );
            }
            (on, _) => panic!("run_subject_units chose the other form (switch {on})"),
        }
    }
}

// No NCBI counterpart: LOSAT_X_SEEDUNMAP's rule at its boundaries: a query-group span or a
// chunk of 65,535 and 65,536 letters, a long side with a short one (the 200 short queries
// against AP027133, and the reverse), no overflow at the largest sizes; never with the switch
// off, and only on Linux with glibc. Two minimum sides make a table of at least 512 KiB, above
// `X_MAPPED_MIN`. The mapped capacities pass glibc's largest dynamic mmap threshold for both
// element types.
#[test]
fn seed_unmap_rule_at_its_boundaries() {
    let glibc = cfg!(all(target_os = "linux", target_env = "gnu"));
    let side = X_MAPPED_MIN_SIDE;
    assert_eq!(side, 65_536);
    assert_eq!(std::mem::size_of::<DiagStruct>(), 4);
    assert!(2 * side * std::mem::size_of::<DiagStruct>() >= X_MAPPED_MIN);
    for (span, chunk_len, holds) in [
        (side - 1, side, false),
        (side, side - 1, false),
        (side, side, true),
        (side - 1, side - 1, false),
        (1 << 22, side, true),
        (side, 1 << 22, true),
        (4_000, 202_000, false),
        (220_000, 1_000, false),
        (220_000, 202_000, true),
        (0, 0, false),
        (usize::MAX, usize::MAX, true),
    ] {
        let at = format!("span {span} chunk {chunk_len}");
        assert_eq!(x_unit_mapped(true, span, chunk_len), glibc && holds, "{at}");
        assert!(!x_unit_mapped(false, span, chunk_len), "{at}: switch off");
    }
    assert_eq!(X_MAPPED_BYTES, (32 << 20) + 4096);
    assert!(x_mapped_len::<DiagStruct>(0) * std::mem::size_of::<DiagStruct>() > X_MAPPED_BYTES);
    assert!(x_mapped_len::<XUnitHsp>(0) * std::mem::size_of::<XUnitHsp>() > X_MAPPED_BYTES);
    assert_eq!(x_mapped_len::<DiagStruct>(1 << 24), 1 << 24);
    assert_eq!(x_mapped_len::<XUnitHsp>(1 << 21), 1 << 21);
}

// No NCBI counterpart: LOSAT_X_SEEDUNMAP's THP gate on the texts of
// /sys/kernel/mm/transparent_hugepage/enabled and /proc/self/status: huge pages round the
// mapped blocks up with "[always]", an unreadable or unrecognised file, unless the process has
// them disabled (`THP_enabled: 0`); not with "[madvise]" or "[never]".
#[test]
fn seed_unmap_thp_gate_reads_the_two_files() {
    let status_1 = "Name:\tLOSAT\nUmask:\t0022\nTHP_enabled:\t1\nThreads:\t8\n";
    let status_0 = "Name:\tLOSAT\nUmask:\t0022\nTHP_enabled:\t0\nThreads:\t8\n";
    let status_without = "Name:\tLOSAT\nUmask:\t0022\nThreads:\t8\n";
    for (enabled, status, rounds_up) in [
        (Some("always [madvise] never\n"), Some(status_1), false),
        (Some("[always] madvise never\n"), Some(status_1), true),
        (Some("always madvise [never]\n"), Some(status_1), false),
        (None, Some(status_1), true),
        (Some("[always] madvise never\n"), Some(status_0), false),
        (None, Some(status_0), false),
        (Some("[always] madvise never\n"), Some(status_without), true),
        (
            Some("always [madvise] never\n"),
            Some(status_without),
            false,
        ),
        (Some("[always] madvise never\n"), None, true),
        (Some("always madvise [never]\n"), None, false),
        (None, None, true),
        (Some(""), None, true),
        (Some("always madvise never\n"), None, true),
        (Some("[madvise"), None, true),
        (
            Some("[always] madvise never\n"),
            Some("THP_enabled:\t01\n"),
            true,
        ),
    ] {
        assert_eq!(
            x_thp_rounds_up(enabled, status),
            rounds_up,
            "{enabled:?} {status:?}"
        );
    }
}

/// One LOSAT_X_SEEDUNMAP case: one query context of `W` blocks of 100 letters around `q_fill`
/// letters of `D`, and two subject frames of `s_len` letters in blocks of 20 `W` and 380 `P`
/// (one chunk each). Every query position has lookup entries (the context's span is about
/// `q_fill + 200`), the subject hits only the `W` blocks, and each block saves about 220 init
/// HSPs, so long sides stay cheap and a unit's copy passes `X_MAPPED_MIN`.
fn unmap_case(q_fill: usize, s_len: usize) -> CompactCase {
    let query = format!("{0}{1}{0}", "W".repeat(100), "D".repeat(q_fill));
    let frames = vec![frame(1, &aa(&query))];
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (lookup, contexts) = build_ncbi_lookup(&[frames], 11, &params, true);
    let block = format!("{}{}", "W".repeat(20), "P".repeat(380));
    let subject: String = block.chars().cycle().take(s_len).collect();
    let shifted: String = block.chars().cycle().skip(200).take(s_len).collect();
    CompactCase {
        pp: XPairPar {
            groups: build_groups(&lookup, &contexts, X_MAX_GROUPS),
            stub: false,
            sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
            sc: Blosum62Score,
            sliced: None,
            mask_rule: None,
            extents: Vec::new(),
            ctx_of: x_ctx_table(&lookup, true),
        },
        cutoff_scores: vec![12; contexts.len()],
        x_dropoff_per_context: vec![16; contexts.len()],
        contexts,
        s_frames: vec![frame(2, &aa(&subject)), frame(-1, &aa(&shifted))],
        failure: crate::algorithm::tblastx::subjects::XReadFailure::new(),
    }
}

/// The units of `case`, run in order on a scratch with LOSAT_X_SEEDUNMAP and on one without:
/// the same records, the same copies (with the mapped capacity where a large copy of a unit of
/// the rule gets it, on Linux with glibc), the same lists of `run_unit` and the same merge.
/// Returns the units where the rule holds and the copies given the mapped capacity.
fn unmap_check(case: &CompactCase) -> (usize, usize) {
    let inp = case.inp();
    let (frame_chunks, units) = list_units(&case.pp, &inp, MAX_DBSEQ_LEN);
    assert!(!units.is_empty());
    let mut today = XScratch::new(false);
    let mut mapped = XScratch::new(true);
    let (mut held, mut copies) = (0usize, 0usize);
    let mut today_records = Vec::new();
    let mut mapped_records = Vec::new();
    for unit in &units {
        let at = format!("unit f {} k {} g {}", unit.f, unit.k, unit.g);
        let span = x_unit_span(&case.pp.groups[unit.g]);
        let rule = x_unit_mapped(true, span, unit.chunk.length);
        held += usize::from(rule);
        let exact = case.pp.run_unit_compact::<false>(unit, &inp, &mut today);
        let copy = case.pp.run_unit_compact::<false>(unit, &inp, &mut mapped);
        assert_eq!(copy, exact, "{at}");
        assert_eq!(
            exact.capacity(),
            exact.len(),
            "{at}: spare room without the switch"
        );
        let large = rule && copy.len() * std::mem::size_of::<XUnitHsp>() >= X_MAPPED_MIN;
        let capacity = if large {
            x_mapped_len::<XUnitHsp>(copy.len())
        } else {
            copy.len()
        };
        assert_eq!(copy.capacity(), capacity, "{at}");
        copies += usize::from(large);
        // Today's lists on the scratch with the mapped table, and the records expand to them.
        let list = case.pp.run_unit::<false>(unit, &inp, &mut today);
        assert_eq!(
            case.pp.run_unit::<false>(unit, &inp, &mut mapped),
            list,
            "{at}"
        );
        let mut expanded = Vec::new();
        x_expand_unit(&copy, unit.f, &inp, &mut expanded);
        assert_eq!(expanded, list, "{at}");
        // LOSAT_X_TBXKMERGE: the task's sort of the copy, with the scratch it was copied from as
        // the buffer (the mapped capacity where the rule holds), is `slice::sort_by` and
        // allocates nothing.
        let buffer = mapped.records.capacity();
        let mut in_task = copy.clone();
        in_task.x_sort_in(&mut mapped);
        let mut sorted = copy.clone();
        sorted.x_sort();
        assert_eq!(in_task, sorted, "{at}");
        assert_eq!(
            mapped.records.capacity(),
            buffer,
            "{at}: the sort allocated"
        );
        today_records.push(exact);
        mapped_records.push(copy);
    }
    if held > 0 {
        assert!(mapped.diag.capacity() >= x_mapped_len::<DiagStruct>(0));
        assert!(mapped.records.capacity() >= x_mapped_len::<XUnitHsp>(0));
    } else {
        assert!(mapped.diag.capacity() < x_mapped_len::<DiagStruct>(0));
        assert!(mapped.records.capacity() < x_mapped_len::<XUnitHsp>(0));
    }
    assert!(today.diag.capacity() < x_mapped_len::<DiagStruct>(0));
    assert!(today.records.capacity() < x_mapped_len::<XUnitHsp>(0));
    let today_hits = XPairPar::<Blosum62Score>::merge_subject_listed(
        XSubjectUnits {
            frame_chunks: frame_chunks.clone(),
            units: units.clone(),
            outputs: XUnitOutputs::Compact(today_records),
        },
        &inp,
    );
    let mapped_hits = XPairPar::<Blosum62Score>::merge_subject_listed(
        XSubjectUnits {
            frame_chunks,
            units,
            outputs: XUnitOutputs::Compact(mapped_records),
        },
        &inp,
    );
    shadow_compare(inp.s_idx, &mapped_hits, &today_hits);
    assert!(!today_hits.is_empty());
    (held, copies)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
// LOSAT_X_SEEDUNMAP on units where the rule holds (a context of about 66,200 letters against
// subject frames of 66,000: both sides long; every unit's copy passes 128 KiB), where the
// query side is short (a context of about 1,200 letters; the copies pass 128 KiB too) and where
// the chunk is short (2,000 letters): the same records, lists and merge, and mapped capacities
// only where the rule holds.
#[test]
fn seed_unmap_units_with_the_mapped_capacities_are_todays() {
    let glibc = cfg!(all(target_os = "linux", target_env = "gnu"));
    for (q_fill, s_len, rule) in [
        (66_000, 66_000, true),
        (1_000, 66_000, false),
        (66_000, 2_000, false),
    ] {
        let at = format!("query fill {q_fill} subjects of {s_len}");
        let case = unmap_case(q_fill, s_len);
        let (held, copies) = unmap_check(&case);
        let expected = if glibc && rule { 2 } else { 0 };
        assert_eq!((held, copies), (expected, expected), "{at}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1043
// On search pools of 1, 2 and 8 threads (at one thread the serial scratch, which never maps;
// at 2 and 8 one scratch per rayon job, with the switch as it acts in this process): the
// records of units where the rule holds are those of a scratch without it.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
#[test]
fn seed_unmap_units_on_search_pools_are_todays() {
    let case = unmap_case(66_000, 66_000);
    let inp = case.inp();
    let (_, units) = list_units(&case.pp, &inp, MAX_DBSEQ_LEN);
    let mut today = XScratch::new(false);
    let expected: Vec<Vec<XUnitHsp>> = units
        .iter()
        .map(|unit| case.pp.run_unit_compact::<false>(unit, &inp, &mut today))
        .collect();
    assert!(expected.iter().all(|records| !records.is_empty()));
    for threads in [1usize, 2, 8] {
        crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
            assert_eq!(pool.enabled(), threads > 1);
            let got = case.pp.run_units(
                &units,
                &inp,
                pool,
                XPairPar::run_unit_compact::<false>,
                false,
            );
            assert_eq!(got, expected, "threads {threads}");
            Ok(())
        })
        .expect("search pool");
    }
}

// LOSAT_X_TBXPREFETCH's rule at its boundaries: off, never; on, only where one of the tables
// has at least the switch's number of KiB (the built-in 1,536 and the smallest value, 2).
#[test]
fn tbx_prefetch_rule_at_its_boundaries() {
    let built_in = X_PREFETCH_MIN_KIB * 1024;
    assert!(!x_prefetch_units(None, [usize::MAX].into_iter()));
    assert!(!x_prefetch_units(
        Some(X_PREFETCH_MIN_KIB),
        std::iter::empty()
    ));
    assert!(!x_prefetch_units(
        Some(X_PREFETCH_MIN_KIB),
        [0, built_in - 1].into_iter()
    ));
    assert!(x_prefetch_units(
        Some(X_PREFETCH_MIN_KIB),
        [0, built_in, 1].into_iter()
    ));
    assert!(x_prefetch_units(Some(2), [2048].into_iter()));
    assert!(!x_prefetch_units(Some(2), [2047].into_iter()));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-516
// The unit functions with the prefetch of LOSAT_X_TBXPREFETCH (`PF`) on search pools of 1, 2
// and 8 threads, in both list forms: the random cases of the compact tests (one and several
// chunks per frame) and the large case of the unmap tests (batches far longer than
// X_PREFETCH_AHEAD) give the lists of the functions without it.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
#[test]
fn tbx_prefetch_units_are_todays() {
    let mut rng = Rng(0x5eed_c0de_0000_0005);
    let mut cases = vec![(unmap_case(66_000, 66_000), MAX_DBSEQ_LEN)];
    for case_no in 0..12 {
        let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
        cases.push(if case_no % 2 == 0 {
            (
                CompactCase::random(&mut rng, &alphabet, 3, 400),
                MAX_DBSEQ_LEN,
            )
        } else {
            (CompactCase::random(&mut rng, &alphabet, 200, 400), 160)
        });
    }
    for threads in [1usize, 2, 8] {
        crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
            let mut records = 0usize;
            for (case_no, (case, max_dbseq_len)) in cases.iter().enumerate() {
                let inp = case.inp();
                let (_, units) = list_units(&case.pp, &inp, *max_dbseq_len);
                let at = format!("threads {threads} case {case_no}");
                let compact = case.pp.run_units(
                    &units,
                    &inp,
                    pool,
                    XPairPar::run_unit_compact::<false>,
                    false,
                );
                records += compact.iter().map(Vec::len).sum::<usize>();
                assert_eq!(
                    case.pp.run_units(
                        &units,
                        &inp,
                        pool,
                        XPairPar::run_unit_compact::<true>,
                        false
                    ),
                    compact,
                    "{at}"
                );
                assert_eq!(
                    case.pp
                        .run_units(&units, &inp, pool, XPairPar::run_unit::<true>, false),
                    case.pp
                        .run_units(&units, &inp, pool, XPairPar::run_unit::<false>, false),
                    "{at}"
                );
            }
            assert!(records > 500, "threads {threads}: records {records}");
            Ok(())
        })
        .expect("search pool");
    }
}

/// What `units_equal_the_masked_reference` saw.
#[derive(Debug)]
struct XMaskedSeen {
    units: usize,
    masked_units: usize,
    saved: usize,
    merged: usize,
    short_frames: usize,
    settings: usize,
    /// Units on masked tables, by matrix, without and with the batched extensions.
    combos: [[usize; 2]; 4],
    /// How often the reference took each branch of the two-hit test.
    outcomes: [usize; 7],
}

// No NCBI counterpart: `-word_size` 2 (SC10b PR-2) and 4 (PR-3). The unit path of a search of
// `word`-letter words (the other matrices' instance, `MatrixScore`, which BLOSUM62 takes at
// those word lengths too) against the reference statements on one masked table
// (`masked_reference`: the table walks of `extend_hit_two_hit_with`), in every instance a unit
// runs: the scan loop and the branch-free one, with and without the prefetch, the sentinel
// walks of 2-letter words (`SENT` = 2, `TbxScore::sent_two_hit_of::<2>`; none at 4 letters)
// and none, the batched extensions, and the masked tables of LOSAT_X_TBXUNITMASK; then the
// subject's merged list on 1, 4 and 8 threads. Query contexts with SEG masks, empty, or of
// one or two letters; subject frames now and then of 0 to 3 letters (shorter than a word, one
// word); X-drops above 127 (the walks refuse them).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1104-1108
/// The comparison above for words of `word` letters: `cases` cases from the generators `seeds`
/// (data, knobs) with the lookup `thresholds`; with `wide`, the units scan the tables of the
/// joint builder (`XPairLutQuery::x_build_tables_wide`, on a pool of 1, 4 or 8 threads by case)
/// in place of `build_groups`' filter of the chain builder's table.
fn units_equal_the_masked_reference(
    word: usize,
    seeds: (u64, u64),
    cases: usize,
    thresholds: [i32; 4],
    wide: bool,
) -> XMaskedSeen {
    use crate::algorithm::tblastx::scoring::MatrixScore;
    let mut rng = Rng(seeds.0);
    let alphabets = [
        aa("AW"),
        aa("ACW"),
        aa("ACDEKW"),
        aa("ACDEFGHIKLMNPQRSTVWY"),
    ];
    let tables: Vec<MatrixScore> = [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Pam30,
        ScoringMatrix::Blosum45,
        ScoringMatrix::Identity,
    ]
    .into_iter()
    .map(|matrix| MatrixScore::with_word(matrix, word))
    .collect();
    // LOSAT_X_TBXUNGAPSENT is on (the batches' contexts have their sentinels).
    let on = x_tbx_sent_contexts(std::iter::empty());
    let (mut units, mut sent_units, mut masked_units) = (0usize, 0usize, 0usize);
    let (mut saved, mut short_frames, mut merged) = (0usize, 0usize, 0usize);
    let mut outcomes = [0usize; 7];
    let serial = crate::utils::threading::SearchPool::x_serial();
    // The sentinel walks take words of 2 and 3 letters (`prepare`).
    let sent_expected = on && matches!(word, 2 | 3);
    // Every knob (matrix, alphabet, window, threshold, masked unit tables) is its own draw
    // from a second generator, so the 300 cases meet in combination (the data generator's
    // stream is left as it was). `combos[matrix][batched]` counts the units on masked tables.
    let mut knobs = Rng(seeds.1);
    let mut combos = [[0usize; 2]; 4];
    let mut settings = std::collections::BTreeSet::new();
    for case_no in 0..cases {
        let matrix_no = knobs.below(tables.len() as u64) as usize;
        let alphabet_no = knobs.below(alphabets.len() as u64) as usize;
        let window_no = knobs.below(4) as usize;
        let threshold_no = knobs.below(4) as usize;
        let masked = knobs.below(2) == 1;
        settings.insert((matrix_no, alphabet_no, window_no, threshold_no, masked));
        let sc = tables[matrix_no];
        assert_eq!(sc.word(), word);
        let alphabet = &alphabets[alphabet_no];
        // The batched extensions only where their rule allows them (`x_tbx_extbatch::active`:
        // `window >= wordsize + BATCH`, here at its bound too); a window of one letter more
        // than a word without them.
        let w = word as i32;
        let window = [40, w + x_tbx_extbatch::BATCH as i32, 64, w + 1][window_no];
        let batched = window >= w + x_tbx_extbatch::BATCH as i32;
        let threshold = thresholds[threshold_no];
        let queries: Vec<Vec<QueryFrame>> = (0..1 + usize::from(case_no % 5 == 4))
            .map(|_| {
                (0..1 + rng.below(6) as usize)
                    .map(|c| random_masked_frame(&mut rng, FRAMES[c], 90, alphabet))
                    .collect()
            })
            .collect();
        let all: Vec<&QueryFrame> = queries.iter().flatten().collect();
        let mut s_frames = Vec::new();
        for &f in &FRAMES[..1 + rng.below(3) as usize] {
            let s_len = if rng.below(6) == 0 {
                short_frames += 1;
                rng.below(4) as usize
            } else {
                60 + rng.below(1200) as usize
            };
            let mut residues = Vec::with_capacity(s_len);
            while residues.len() < s_len {
                let q = all[rng.below(all.len() as u64) as usize];
                if rng.below(3) != 0 && q.aa_len > 0 {
                    let from = rng.below(q.aa_len as u64) as usize;
                    let take = (q.aa_len - from).min(s_len - residues.len());
                    residues.extend_from_slice(&q.aa_seq[1 + from..1 + from + take]);
                } else {
                    let take = (2 + rng.below(30) as usize).min(s_len - residues.len());
                    residues.extend(random_residues(&mut rng, take, alphabet));
                }
            }
            s_frames.push(frame(f, &residues));
        }
        let params = lookup_protein_params_ungapped(sc.matrix());
        let (lookup, contexts) = crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile(
            &queries,
            threshold,
            &params,
            true,
            sc.matrix(),
            word,
        );
        assert_eq!(lookup.word_length, word as i32);
        let query_length = query_length_of(&contexts);
        let size = diag_array_size_for(query_length, window);
        let cutoff_scores: Vec<i32> = (0..contexts.len())
            .map(|i| 8 + (i as i32 % 4) * 3)
            .collect();
        let x_dropoff_per_context: Vec<i32> = (0..contexts.len())
            .map(|i| [7, 16, 30, 200][(i + case_no) % 4])
            .collect();
        let mut diag_array = vec![DiagStruct::default(); size as usize];
        let mut diag_offset = window;
        let (reference, reference_outcomes, _) = masked_reference(
            sc,
            &lookup,
            &contexts,
            &s_frames,
            &cutoff_scores,
            &x_dropoff_per_context,
            window,
            &mut diag_array,
            &mut diag_offset,
        );
        for outcome in reference_outcomes {
            outcomes[outcome as usize] += 1;
        }
        let inp = XSubject {
            s_frames_preliminary: &s_frames,
            s_frames: &s_frames,
            contexts: &contexts,
            cutoff_scores: &cutoff_scores,
            x_dropoff_per_context: &x_dropoff_per_context,
            window,
            wordsize: word as i32,
            s_idx: 0,
            s_len: 999,
            failure: &crate::algorithm::tblastx::subjects::XReadFailure::new(),
            passed: None,
            subject_besthit: false,
        };
        // Masked unit tables in about half of the cases (ratio 1: every unit whose unmasked
        // table is at least the masked one's size).
        // The joint builder's tables on a pool of 1, 4 or 8 threads (by case), or the filter.
        let groups = if wide {
            let q = XPairLutQuery::x_prepare_for(&queries, threshold, true, sc.matrix(), word)
                .expect("the joint builder covers the word length");
            let locations = q
                .x_group_locations(&group_starts(&contexts, X_MAX_GROUPS))
                .expect("every location inside one group");
            let build = |pool: &crate::utils::threading::SearchPool<'_>| {
                q.x_build_tables_wide(&locations, pool, true)
                    .into_iter()
                    .map(XGroup::from_direct)
                    .collect::<Vec<_>>()
            };
            match case_no % 3 {
                0 => build(&serial),
                threads => {
                    crate::utils::threading::with_search_pool(4 * threads, "tblastx", |pool| {
                        Ok(build(pool))
                    })
                    .expect("search pool")
                }
            }
        } else {
            build_groups(&lookup, &contexts, X_MAX_GROUPS)
        };
        let pp = XPairPar {
            groups,
            stub: false,
            // As `prepare` decides it.
            sent: matches!(sc.word(), 2 | 3)
                && sc.x_sent_walks()
                && x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
            sc,
            sliced: None,
            mask_rule: masked.then_some(XMaskRule {
                max_cells: usize::MAX,
                ratio: 1,
                band_cells: None,
            }),
            extents: if masked {
                x_group_extents(&contexts, X_MAX_GROUPS)
            } else {
                Vec::new()
            },
            ctx_of: x_ctx_table(&lookup, true),
        };
        assert_eq!(pp.sent, sent_expected, "case {case_no}");
        assert_eq!(
            unit_lists(&pp, &inp),
            reference,
            "case {case_no}: unit lists differ from the masked reference"
        );
        if batched {
            assert_eq!(
                unit_lists_batched(&pp, &inp),
                reference,
                "case {case_no}: batched unit lists differ from the masked reference"
            );
        }
        saved += reference.iter().map(Vec::len).sum::<usize>();
        // Every instance of a unit, against the scan loop without the walks.
        let mut scratch = XScratch::new(false);
        for (f, s_frame) in s_frames.iter().enumerate() {
            let chunk = SubjectChunk {
                offset: 0,
                length: s_frame.aa_len,
                overlap: 0,
            };
            for g in 0..pp.groups.len() {
                if pp.groups[g].lookup.longest_chain <= 0 {
                    continue;
                }
                let unit = pp.x_unit(&inp, f, 0, chunk, g);
                let sent = pp.x_unit_sent(&unit, &inp);
                assert_eq!(sent, pp.sent);
                units += 1;
                sent_units += usize::from(sent);
                let is_masked = pp.x_unit_mask(&unit, window).is_some();
                masked_units += usize::from(is_masked);
                combos[matrix_no][usize::from(batched)] += usize::from(is_masked);
                let mut want: Vec<InitHSP> = Vec::new();
                pp.scan_unit_hits::<_, XNoProbe, false, false, 0>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut want,
                    &mut XNoProbe,
                );
                let mut lists: Vec<(&str, Vec<InitHSP>)> = Vec::new();
                let mut run = |name, got: Vec<InitHSP>| lists.push((name, got));
                let mut got = Vec::new();
                pp.scan_unit_hits::<_, XNoProbe, true, true, 0>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut got,
                    &mut XNoProbe,
                );
                run("bf pf", got);
                if batched {
                    let mut got = Vec::new();
                    pp.scan_unit_batched_hits::<_, true, true, 0>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        x_tbx_extbatch::BATCH,
                    );
                    run("batched bf pf", got);
                }
                if sent {
                    for (name, pf, bf) in [
                        ("sent", false, false),
                        ("bf sent", false, true),
                        ("bf pf sent", true, true),
                    ] {
                        let mut got = Vec::new();
                        match (pf, bf) {
                            (false, false) => pp.scan_unit_hits::<_, XNoProbe, false, false, 2>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                                &mut XNoProbe,
                            ),
                            (false, true) => pp.scan_unit_hits::<_, XNoProbe, false, true, 2>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                                &mut XNoProbe,
                            ),
                            _ => pp.scan_unit_hits::<_, XNoProbe, true, true, 2>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                                &mut XNoProbe,
                            ),
                        }
                        run(name, got);
                    }
                    for len in [x_tbx_extbatch::BATCH, 1].into_iter().filter(|_| batched) {
                        let mut got = Vec::new();
                        pp.scan_unit_batched_hits::<_, true, true, 2>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            len,
                        );
                        run("batched bf pf sent", got);
                        let mut got = Vec::new();
                        pp.scan_unit_batched_hits::<_, false, false, 2>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            len,
                        );
                        run("batched sent", got);
                    }
                }
                for (name, got) in &lists {
                    assert_eq!(got, &want, "case {case_no} frame {f} group {g}: {name}");
                }
            }
        }
        // The subject's merged list on 1, 4 and 8 threads (`run_subject_units`, whose units
        // take `run_unit` or its batched form by the rules, and `merge_subject`).
        let one = XPairPar::<MatrixScore>::merge_subject(pp.run_subject_units(&inp, &serial), &inp);
        merged += one.len();
        // Against the reference's own lists, converted and merged by the reference's
        // statements (one chunk per frame here: `get_ungapped_hsp_list` per frame list, the
        // frame lists appended and sorted), not by anything of `run_subject_units`.
        let mut expect: Vec<UngappedHit> = Vec::new();
        for list in &reference {
            if !list.is_empty() {
                expect.extend(get_ungapped_hsp_list(
                    list.clone(),
                    inp.contexts,
                    inp.s_frames,
                ));
                if !ungapped_hits_is_sorted_by_score_ncbi(&expect) {
                    sort_ungapped_hits_by_score_ncbi(&mut expect);
                }
            }
        }
        assert_eq!(
            one.len(),
            expect.len(),
            "case {case_no}: merged list length"
        );
        for (k, (a, b)) in one.iter().zip(&expect).enumerate() {
            assert!(
                x_same_hit(a, b),
                "case {case_no}: merged hit {k} (1 thread)"
            );
        }
        for threads in [4usize, 8] {
            crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
                let got =
                    XPairPar::<MatrixScore>::merge_subject(pp.run_subject_units(&inp, pool), &inp);
                assert_eq!(got.len(), expect.len(), "case {case_no} threads {threads}");
                for (k, (a, b)) in got.iter().zip(&expect).enumerate() {
                    assert!(x_same_hit(a, b), "case {case_no} threads {threads} hit {k}");
                }
                Ok(())
            })
            .expect("search pool");
        }
    }
    assert_eq!(sent_units, if sent_expected { units } else { 0 });
    XMaskedSeen {
        units,
        masked_units,
        saved,
        merged,
        short_frames,
        settings: settings.len(),
        combos,
        outcomes,
    }
}

// No NCBI counterpart: `-word_size 2` (SC10b PR-2). The unit path of a search of 2-letter
// words against the reference statements on one masked table (`units_equal_the_masked_reference`).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1104-1108
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn units_of_2_letter_words_equal_the_masked_reference() {
    let seen = units_equal_the_masked_reference(
        2,
        (0x2b17_5eed_0000_0002, 0x2b17_6b0b_0000_0002),
        300,
        [11, 13, 0, 8],
        false,
    );
    assert_masked_seen(&seen, (1_000, 350, 100_000, 100_000, 100, 200));
}

// No NCBI counterpart: `-word_size 4` (SC10b PR-3). The unit path of a search of 4-letter
// words, on the tables of the joint builder (`x_build_tables_wide`, on pools of 1, 4 and 8
// threads), against the reference statements on one masked table over the chain builder's
// full table (`units_equal_the_masked_reference`): the scan loop and the branch-free one, with
// and without the prefetch, the batched extensions, the masked tables of LOSAT_X_TBXUNITMASK
// (no sentinel walks at 4 letters), then the merged list on 1, 4 and 8 threads.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1104-1108
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn units_of_4_letter_words_equal_the_masked_reference() {
    let seen = units_equal_the_masked_reference(
        4,
        (0x4b17_5eed_0000_0004, 0x4b17_6b0b_0000_0004),
        200,
        [11, 13, 0, 15],
        true,
    );
    assert_masked_seen(&seen, (700, 200, 20_000, 20_000, 60, 150));
}

/// The floors of `units_equal_the_masked_reference`'s counts: (units, masked units, saved,
/// merged, short frames, distinct settings); every matrix met masked tables with and without the
/// batched extensions, and the reference took every branch of the two-hit test.
fn assert_masked_seen(seen: &XMaskedSeen, floors: (usize, usize, usize, usize, usize, usize)) {
    let (units, masked_units, saved, merged, short_frames, settings) = floors;
    assert!(seen.units > units, "{seen:?}");
    assert!(seen.masked_units > masked_units, "{seen:?}");
    assert!(seen.saved > saved, "{seen:?}");
    assert!(seen.merged > merged, "{seen:?}");
    assert!(seen.short_frames > short_frames, "{seen:?}");
    assert!(seen.settings > settings, "{seen:?}");
    for (m, row) in seen.combos.iter().enumerate() {
        for (b, n) in row.iter().enumerate() {
            assert!(*n > 0, "no masked unit: matrix {m}, batched {b}");
        }
    }
    for (outcome, n) in seen.outcomes.iter().enumerate() {
        assert!(*n > 0, "outcome {outcome} never taken");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456,461-465
// `-word_size 2` (SC10b PR-2): the direct groups of 2-letter words (10-bit backbone indices)
// are `build_groups` of the reference builder's table of that word length, for one query and
// several (groups of several contexts), SEG masks and contexts without words, thresholds 0
// (no neighbours), 8, 11, 13 and INT_MIN (every word), several matrices, and on a search pool
// (the memo of neighbour lists keyed by the word length). Words of 4 letters (20-bit indices)
// take the joint builder (`wide_groups_of_4_letter_words_equal_build_groups`).
#[test]
fn direct_groups_of_2_letter_words_equal_build_groups() {
    use crate::algorithm::tblastx::lookup::x_pair_lut::x_direct_covers;
    assert!(x_direct_covers(2) && x_direct_covers(3) && !x_direct_covers(4));
    let mut rng = Rng(0x2b17_d1ec_0000_0002);
    let alphabets = [aa("ACDEKW"), aa("AW"), aa("ACDEFGHIKLMNPQRSTVWY")];
    let matrices = [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Pam30,
        ScoringMatrix::Blosum45,
        ScoringMatrix::Pam250,
    ];
    let serial = crate::utils::threading::SearchPool::x_serial();
    let (mut with_entries, mut empty_groups, mut six_groups) = (0usize, 0usize, 0usize);
    for case_no in 0..48 {
        let alphabet = &alphabets[case_no % alphabets.len()];
        let matrix = matrices[case_no % matrices.len()];
        let queries: Vec<Vec<QueryFrame>> = (0..if case_no % 3 == 2 { 3 } else { 1 })
            .map(|_| {
                let n_ctx = 1 + rng.below(6) as usize;
                (0..n_ctx)
                    .map(|c| random_masked_frame(&mut rng, FRAMES[c], 200, alphabet))
                    .collect()
            })
            .collect();
        for threshold in [0, 8, 11, 13] {
            let (groups, filled) = lut_check_for(&queries, threshold, &serial, matrix, 2);
            with_entries += filled;
            empty_groups += groups - filled;
            six_groups += usize::from(groups == X_MAX_GROUPS);
        }
    }
    assert!(with_entries > 300, "groups with entries {with_entries}");
    assert!(empty_groups > 50, "empty groups {empty_groups}");
    assert!(six_groups > 10, "six groups {six_groups}");
    // INT_MIN: every 2-letter word is a neighbour of every query word.
    let alphabet = aa("ACDEFGHIKLMNPQRSTVWY");
    let frames: Vec<QueryFrame> = FRAMES
        .iter()
        .map(|&f| random_masked_frame(&mut rng, f, 30, &alphabet))
        .collect();
    assert_eq!(
        lut_check_for(
            std::slice::from_ref(&frames),
            i32::MIN,
            &serial,
            ScoringMatrix::Blosum62,
            2
        )
        .0,
        6
    );
    // On a search pool, twice (the second round reads the memoised lists), and with a build of
    // 3-letter words between (another key of the memo).
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    crate::utils::threading::with_search_pool(4, "tblastx", |pool| {
        for _round in 0..2 {
            for word in [2, 3, 2] {
                let (groups, _) = lut_check_for(
                    std::slice::from_ref(&frames),
                    13,
                    pool,
                    ScoringMatrix::Pam30,
                    word,
                );
                assert_eq!(groups, 6);
            }
        }
        Ok(())
    })
    .expect("search pool");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456,461-465
// `-word_size 4` (SC10b PR-3): the groups of the joint builder (`lut_wide_groups`, 20-bit
// backbone indices) are `build_groups` of the chain builder's table of that word length, with the
// stub, for one query and several (groups of several contexts), SEG masks and contexts without
// words, thresholds 0, 11, 13 and INT_MIN, several matrices, on one thread and on pools of 4
// and 8 threads.
#[test]
fn wide_groups_of_4_letter_words_equal_build_groups() {
    let mut rng = Rng(0x4b17_d1ec_0000_0004);
    let alphabets = [aa("ACDEKW"), aa("AW"), aa("ACDEFGHIKLMNPQRSTVWY")];
    let matrices = [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Pam30,
        ScoringMatrix::Blosum62,
        ScoringMatrix::Pam250,
    ];
    let cases: Vec<(Vec<Vec<QueryFrame>>, ScoringMatrix)> = (0..24)
        .map(|case_no| {
            let alphabet = &alphabets[case_no % alphabets.len()];
            let queries: Vec<Vec<QueryFrame>> = (0..if case_no % 3 == 2 { 3 } else { 1 })
                .map(|_| {
                    let n_ctx = 1 + rng.below(6) as usize;
                    (0..n_ctx)
                        .map(|c| random_masked_frame(&mut rng, FRAMES[c], 120, alphabet))
                        .collect()
                })
                .collect();
            (queries, matrices[case_no % matrices.len()])
        })
        .collect();
    let alphabet = aa("ACDEFGHIKLMNPQRSTVWY");
    let tiny: Vec<QueryFrame> = FRAMES
        .iter()
        .map(|&f| random_masked_frame(&mut rng, f, 5, &alphabet))
        .collect();
    let check = |pool: &crate::utils::threading::SearchPool<'_>| {
        let (mut with_entries, mut empty_groups, mut six_groups) = (0usize, 0usize, 0usize);
        for (queries, matrix) in &cases {
            for threshold in [0, 11, 13] {
                let (groups, filled) = lut_check_for(queries, threshold, pool, *matrix, 4);
                with_entries += filled;
                empty_groups += groups - filled;
                six_groups += usize::from(groups == X_MAX_GROUPS);
            }
        }
        assert!(with_entries > 100, "groups with entries {with_entries}");
        assert!(empty_groups > 10, "empty groups {empty_groups}");
        assert!(six_groups > 5, "six groups {six_groups}");
        // INT_MIN: every 4-letter word is a neighbour of every query word.
        assert_eq!(
            lut_check_for(
                std::slice::from_ref(&tiny),
                i32::MIN,
                pool,
                ScoringMatrix::Blosum62,
                4
            )
            .0,
            6
        );
    };
    check(&crate::utils::threading::SearchPool::x_serial());
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    for threads in [4usize, 8] {
        crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
            check(pool);
            Ok(())
        })
        .expect("search pool");
    }
}

// No NCBI counterpart: LOSAT_X_TBXBAND's memory rule with LOSAT_X_TBXUNITMASK (`x_band_masked`)
// for a batch whose groups are never cut (words of 4 letters, `x_slices_fit`): where
// LOSAT_X_TBXQSLICE's rule would cut the group (200,000 letters at 32 threads), `decide` keeps the
// guard and `decide_uncut` takes the units, as `prepare` keeps the masked tables of such a batch;
// `decide_for` chooses by the word length.
#[test]
fn groups_never_cut_keep_the_masked_rule() {
    let mut rng = Rng(13);
    let (_, mut contexts) = lookup_for(&[frame(1, &random_residues(&mut rng, 216, &aa("ACDEKW")))]);
    contexts.truncate(1);
    contexts[0].aa_len = 200_000;
    let query_length = query_length_of(&contexts);
    let diag_array_size = diag_array_size_for(query_length, 40);
    assert_eq!(diag_array_size, 1 << 18);
    // A subject chunk longer than the table (the first guard fails, `x_band_fits` too).
    let long = [FastaRecord::new("", b"l", &vec![b'A'; 3 * 300_000])];
    let threads = 32;
    let subjects = crate::algorithm::tblastx::subjects::Subjects::Fasta(&long);
    let starts = group_starts(&contexts, x_groups_wanted(threads, subjects));
    // The slicing rule (x86-64 only) cuts this group at 32 threads.
    if !cfg!(target_arch = "x86_64") {
        return;
    }
    assert!(x_tbx_qslice::rule_cuts(&contexts, &starts, threads));
    let gates = XGates {
        contexts: &contexts,
        query_length,
        diag_array_size,
        subjects,
        reference_loop: false,
        chunk_switches: false,
        threads,
        cost_gate: false,
        n1_gate: false,
        threshold: 13,
        band: true,
        unit_mask: true,
    };
    assert_eq!(decide(&gates), XDecision::Guard);
    assert_eq!(decide_uncut(&gates), XDecision::Units);
    assert_eq!(decide_for(3, &gates), XDecision::Guard);
    assert_eq!(decide_for(2, &gates), XDecision::Guard);
    assert_eq!(decide_for(4, &gates), XDecision::Units);
    // Without LOSAT_X_TBXUNITMASK both keep the guard; the gates still come first.
    let unmasked = XGates {
        unit_mask: false,
        ..gates
    };
    assert_eq!(decide_uncut(&unmasked), XDecision::Guard);
    let gated = XGates {
        reference_loop: true,
        ..gates
    };
    assert_eq!(decide_uncut(&gated), XDecision::Gate);
    // Against `decide`'s statements with `x_band_masked` not reading the slicing rule, over
    // gates drawn at random (thread counts, switches, context and subject lengths, tables).
    let reference = |g: &XGates<'_>| -> XDecision {
        if g.reference_loop || g.chunk_switches {
            return XDecision::Gate;
        }
        let (max_context_len, longest_subject_chunk) = x_guard_lengths(g.contexts, g.subjects);
        let first_guard = g.band
            && (x_band_fits(longest_subject_chunk, g.diag_array_size) || g.unit_mask)
            || guard(
                max_context_len,
                longest_subject_chunk,
                g.query_length,
                g.diag_array_size,
            );
        if !first_guard || !last_hit_fits(longest_subject_chunk, g.query_length, g.diag_array_size)
        {
            return XDecision::Guard;
        }
        if g.cost_gate && !units_gain(g) {
            return XDecision::Cost;
        }
        XDecision::Units
    };
    let mut rng = Rng(0x4b17_dec1_0000_0004);
    let mut seen = std::collections::BTreeSet::new();
    for _ in 0..2_000 {
        let mut drawn = contexts.clone();
        drawn[0].aa_len = [100, 5_000, 70_000, 200_000, 1_500_000][rng.below(5) as usize];
        let query_length = query_length_of(&drawn);
        let window = [40, 4_000, 300_000][rng.below(3) as usize];
        let subject = vec![b'A'; [30, 3_000, 900_000][rng.below(3) as usize]];
        let records = [FastaRecord::new("", b"s", &subject)];
        let g = XGates {
            contexts: &drawn,
            query_length,
            diag_array_size: diag_array_size_for(query_length, window),
            subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(&records),
            reference_loop: rng.below(8) == 0,
            chunk_switches: rng.below(8) == 0,
            threads: [1, 4, 8, 32][rng.below(4) as usize],
            cost_gate: rng.below(2) == 0,
            n1_gate: rng.below(2) == 0,
            threshold: [11, 13][rng.below(2) as usize],
            band: rng.below(4) != 0,
            unit_mask: rng.below(4) != 0,
        };
        let got = decide_uncut(&g);
        assert_eq!(got, reference(&g));
        seen.insert(format!("{got:?}"));
    }
    assert_eq!(seen.len(), 4, "every decision met: {seen:?}");
}

// No NCBI counterpart: the slicing rule's table term (`x_tbx_qslice::x_slices_fit`): words of 2
// and 3 letters may be cut, words of 4 letters are not.
#[test]
fn groups_of_4_letter_words_are_not_cut() {
    assert!(x_tbx_qslice::x_slices_fit(2) && x_tbx_qslice::x_slices_fit(3));
    assert!(!x_tbx_qslice::x_slices_fit(4));
}
