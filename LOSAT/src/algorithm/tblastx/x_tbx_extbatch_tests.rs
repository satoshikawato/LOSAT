use super::*;

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

/// Residues from `alphabet`, with stretches copied from `source` (shared diagonals: long
/// extensions, right extensions and later hits on their cells).
fn residues(rng: &mut Rng, len: usize, alphabet: &[u8], source: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    while out.len() < len {
        if !source.is_empty() && rng.below(3) == 0 {
            let at = rng.below(source.len() as u64) as usize;
            let take = (rng.below(60) as usize + 5).min(source.len() - at);
            out.extend_from_slice(&source[at..at + take]);
        } else {
            out.push(alphabet[rng.below(alphabet.len() as u64) as usize]);
        }
    }
    out.truncate(len);
    out
}

// No NCBI counterpart: compares two LOSAT scans of the same units (module doc, items 1-3).
#[test]
fn batched_extensions_equal_scan_unit() {
    let alphabet: Vec<u8> = "ACDEKLMNPQRSTVWY"
        .bytes()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect();
    let narrow: Vec<u8> = "AKLE"
        .bytes()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect();
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let mut rng = Rng(0x5eed_a5a5_0123_4567);
    let mut total_cands = 0u64;
    let mut total_early = 0u64;
    for case in 0..16 {
        let abc = if case % 3 == 2 { &narrow } else { &alphabet };
        let q_len = 200 + rng.below(400) as usize;
        let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
            .iter()
            .map(|&f| frame(f, &residues(&mut rng, q_len, abc, &[])))
            .collect();
        let source: Vec<u8> = frames
            .iter()
            .flat_map(|f| f.aa_seq[1..f.aa_len + 1].to_vec())
            .collect();
        let s_len = 300 + rng.below(500) as usize;
        // Every fourth case scans the query against itself: long walks, walks that reach the
        // ends of both sequences.
        let s_frames: Vec<QueryFrame> = if case % 4 == 3 {
            frames.clone()
        } else {
            [1i8, 2, 3, -1, -2, -3]
                .iter()
                .map(|&f| frame(f, &residues(&mut rng, s_len, abc, &source)))
                .collect()
        };
        let (lookup, contexts) =
            build_ncbi_lookup(std::slice::from_ref(&frames), 11, &params, true);
        let cutoff_scores: Vec<i32> = (0..contexts.len())
            .map(|i| 10 + (i as i32 % 4) * 3)
            .collect();
        // X-drops inside and outside the 1..=127 of the two-hit table walks (their fallback).
        let x_dropoff_per_context: Vec<i32> = (0..contexts.len())
            .map(|i| [7, 16, 30, 200][i % 4])
            .collect();
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
        // Windows the module takes (`active`): from its smallest, wordsize + BATCH.
        let window = [40, 3 + BATCH as i32, 64][case % 3];
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
        let mut scratch = XScratch::new(false);
        for (f, s_frame) in s_frames.iter().enumerate() {
            // Every other case scans a chunk that starts inside the frame (offset > 0).
            let offset = if case % 2 == 1 { s_frame.aa_len / 3 } else { 0 };
            let chunk = SubjectChunk {
                offset,
                length: s_frame.aa_len - offset,
                overlap: 0,
            };
            for (g, group) in pp.groups.iter().enumerate() {
                if group.lookup.longest_chain <= 0 {
                    continue;
                }
                let unit = pp.x_unit(&inp, f, 0, chunk, g);
                let mut reference: Vec<InitHSP> = Vec::new();
                pp.scan_unit::<_, _, false>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut reference,
                    &mut XNoProbe,
                );
                let reference_diag = scratch.diag.clone();
                for len in [1, 2, 3, 8, 16] {
                    let mut batched: Vec<InitHSP> = Vec::new();
                    let (cands, early) = if len == 8 {
                        pp.scan_unit_batched::<_, true>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut batched,
                            len,
                        )
                    } else {
                        pp.scan_unit_batched::<_, false>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut batched,
                            len,
                        )
                    };
                    assert_eq!(
                        batched, reference,
                        "case {case} frame {f} group {g} batch {len}: records"
                    );
                    assert!(
                        scratch.diag == reference_diag,
                        "case {case} frame {f} group {g} batch {len}: table"
                    );
                    if len == 16 {
                        total_cands += cands;
                        total_early += early;
                    }
                }
                let mut compact: Vec<XUnitHsp> = Vec::new();
                let mut compact_ref: Vec<XUnitHsp> = Vec::new();
                pp.scan_unit::<_, _, false>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut compact_ref,
                    &mut XNoProbe,
                );
                pp.scan_unit_batched::<_, false>(&unit, &inp, &mut scratch, &mut compact, 5);
                assert_eq!(
                    compact, compact_ref,
                    "case {case} frame {f} group {g}: compact records"
                );
            }
        }
    }
    // The cases reach both resolution paths: batches and hits on pending cells.
    assert!(total_cands > 1000, "candidates {total_cands}");
    assert!(
        total_early > 10,
        "early resolutions {total_early} of {total_cands}"
    );
}

// No NCBI counterpart: the size rule of the switch (module doc).
#[test]
fn the_rule_takes_the_largest_span() {
    let spans = [100 * 1024, 300 * 1024 - 1, 50];
    assert!(!x_extbatch_units(None, spans.into_iter()));
    assert!(x_extbatch_units(Some(256), spans.into_iter()));
    assert!(x_extbatch_units(Some(299), spans.into_iter()));
    assert!(!x_extbatch_units(Some(300), spans.into_iter()));
    assert!(x_extbatch_units(Some(2), [2048].into_iter()));
    assert!(!x_extbatch_units(Some(2), [2047].into_iter()));
    assert!(!x_extbatch_units(Some(2), std::iter::empty()));
}

// No NCBI counterpart: the batch resolution against `extend_hit_two_hit_fitting` on every candidate,
// with the cell writes and the init HSPs of `x_hit`, some candidates resolved early.
#[test]
fn batch_resolution_equals_extend_hit_two_hit() {
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let mut rng = Rng(0x0dd_ba11_cafe_f00d);
    let mut cands = 0u64;
    for round in 0..400 {
        let q_len = 3 + rng.below(if round % 5 == 0 { 40 } else { 300 }) as usize;
        // Codes 1..=24, a code 0 (BLAST_SCORE_MIN) now and then.
        let pick = |rng: &mut Rng| {
            if rng.below(97) == 0 {
                0u8
            } else {
                1 + rng.below(24) as u8
            }
        };
        let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
            .iter()
            .map(|&f| frame(f, &(0..q_len).map(|_| pick(&mut rng)).collect::<Vec<u8>>()))
            .collect();
        let (lookup, contexts) =
            build_ncbi_lookup(std::slice::from_ref(&frames), 11, &params, true);
        // The subject: pieces of the query frames with substitutions, and random residues.
        let mut subject: Vec<u8> = Vec::new();
        let s_len = 3 + rng.below(400) as usize;
        while subject.len() < s_len {
            if rng.below(2) == 0 {
                let f = &frames[rng.below(6) as usize].aa_seq[1..q_len + 1];
                let at = rng.below(q_len as u64) as usize;
                let take = (1 + rng.below(80) as usize).min(q_len - at);
                let rate = [3u64, 8, 30, 1_000_000][rng.below(4) as usize];
                subject.extend(f[at..at + take].iter().map(|&c| {
                    if rng.below(rate) == 0 {
                        pick(&mut rng)
                    } else {
                        c
                    }
                }));
            } else {
                subject.push(pick(&mut rng));
            }
        }
        subject.truncate(s_len);
        let cutoff_scores: Vec<i32> = (0..contexts.len()).map(|i| 8 + i as i32).collect();
        let x_dropoff_per_context: Vec<i32> = (0..contexts.len()).map(|_| 0).collect();
        let window = 40;
        let failure = crate::algorithm::tblastx::subjects::XReadFailure::new();
        let ctx_of = x_ctx_table(&lookup, true);
        let h = XHitCtx {
            diag_offset: window,
            window,
            wordsize: 3,
            lookup: &lookup,
            ctx_of: &ctx_of,
            contexts: &contexts,
            cutoff_scores: &cutoff_scores,
            x_dropoff_per_context: &x_dropoff_per_context,
            subject: &subject,
            subject_sent: &[],
            s_f_idx: 2,
            s_frame: -1,
            s_idx: 7,
            s_len: 999,
            sc: Blosum62Score,
            query_batch: QueryBatch::new(&contexts, &failure),
        };
        let n = 1 + rng.below(BATCH as u64) as usize;
        let mut b = XPendBatch {
            c: [XCand::default(); BATCH],
            n,
            len: BATCH,
            cands: 0,
            early: 0,
        };
        let mut diag = vec![DiagStruct::default(); n];
        let mut expect_cells = Vec::new();
        let mut expect_hsps: Vec<InitHSP> = Vec::new();
        for i in 0..n {
            let ctx_idx = rng.below(contexts.len() as u64) as usize;
            let ctx = &contexts[ctx_idx];
            let q = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
            // A word that fits (a lookup hit), as `two_hit.rs` reads it unchecked.
            let q_raw = rng.below((q.len() - 2) as u64) as usize;
            let s_right = rng.below((subject.len() - 2) as u64) as usize;
            let diff = 3 + rng.below(37) as usize;
            let s_left = (s_right + 3).saturating_sub(diff);
            let x_drop = [1, 2, 7, 16, 30, 127, 128, 300][rng.below(8) as usize];
            b.c[i] = XCand {
                query_offset: ctx.frame_base as u32 + q_raw as u32,
                s_right: s_right as u32,
                s_left: s_left as u32,
                q_raw: q_raw as u32,
                cell: i as u32,
                ctx_idx: ctx_idx as u32,
                x_drop,
                ..XCand::default()
            };
            diag[i].set_raw_bits(FLAG_BIT | i as u32);
            let (q0, qe, s0, _se, score, right, s_last) = extend_hit_two_hit_fitting(
                q,
                crate::algorithm::tblastx::extension::InBatch::ALONE,
                &subject,
                s_left,
                s_right,
                q_raw,
                x_drop,
                false,
            );
            let mut cell = DiagStruct::default();
            if right {
                cell.set_flag(1);
                cell.set_last_hit(s_last as i32 - 2 + window);
            } else {
                cell.set_last_hit(s_right as i32 + window);
            }
            expect_cells.push(cell);
            if score >= cutoff_scores[ctx_idx] {
                let (q0, qe, s0) = (q0 as i32, qe as i32, s0 as i32);
                expect_hsps.push(InitHSP {
                    q_start_absolute: ctx.frame_base + q0,
                    q_end_absolute: ctx.frame_base + qe,
                    s_start: s0,
                    s_end: s0 + (qe - q0),
                    q_seed_absolute: ctx.frame_base + q_raw as i32,
                    s_seed: s_right as i32,
                    score,
                    ctx_idx,
                    s_f_idx: 2,
                    q_idx: ctx.q_idx,
                    s_idx: 7,
                    q_frame: ctx.frame,
                    s_frame: -1,
                    q_orig_len: ctx.orig_len,
                    s_orig_len: 999,
                });
            }
        }
        // A few candidates resolved already, as a hit on their cell does.
        for i in 0..n {
            if rng.below(10) == 0 {
                // SAFETY: the cell is below n = diag.len().
                unsafe { resolve::<_, 0>(&mut b.c[i], diag.as_mut_ptr(), &h) };
            }
        }
        let mut got: Vec<InitHSP> = Vec::new();
        // SAFETY: every candidate's cell is below n = diag.len().
        unsafe { flush::<_, _, 0>(&mut b, diag.as_mut_ptr(), &h, &mut got) };
        assert_eq!(diag, expect_cells, "round {round}: cells");
        assert_eq!(got, expect_hsps, "round {round}: init HSPs");
        assert_eq!(b.n, 0);
        cands += n as u64;
    }
    assert!(cands > 1000, "candidates {cands}");
}
