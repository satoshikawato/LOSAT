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

/// xorshift64*, for reproducible random values.
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

/// A probe that counts the outcomes of the reference loop's hits.
#[derive(Default)]
struct XCount([u64; 7]);

impl XProbe for XCount {
    type Mark = ();
    fn unit_begin(&mut self, _width: usize) {}
    fn batch(&mut self, _hits: usize) {}
    fn at(&mut self, _cell: usize, _s: u32) {}
    fn hit(&mut self, outcome: XHit, _window: i32) {
        self.0[outcome as usize] += 1;
    }
    fn batch_end(&mut self) {}
    fn ext_begin(&mut self) {}
    fn ext_end<M: TbxScore>(
        &mut self,
        _mark: (),
        _sc: M,
        _q: &[u8],
        _s: &[u8],
        _s_left: usize,
        _s_right: usize,
        _q_right: usize,
        _x_drop: i32,
    ) {
    }
    fn unit_end(&mut self) {}
}

/// What the tests' `rare` does to a cell (the same in both loops): a flag-1 cell is masked
/// below its `last_hit` and reset above it, a candidate sets the flag with `last_hit = s +
/// window` on odd s and stores `s + window` on even s, as `x_hit` stores one of those.
fn rare_apply(cell: &mut DiagStruct, s: u32, window: i32) {
    let next = s.wrapping_add(window as u32) as i32;
    if cell.flag() != 0 {
        if (next as u32) >= cell.last_hit() as u32 {
            cell.set_last_hit(next);
            cell.set_flag(0);
        }
    } else {
        cell.set_last_hit(next);
        cell.set_flag(s & 1);
    }
}

// No NCBI counterpart: the test's per-hit body, `rare_apply` on the hit's cell, and for every
// third subject offset a write to another cell (the next one, cyclically) of the value a
// too-far write stores, `s + window`, as the batched body writes other cells when it resolves
// held-back extensions (x_tbx_extbatch.rs).
/// SAFETY: `table` points to `width` cells and `cell < width`.
unsafe fn rare_and_other(table: *mut DiagStruct, width: usize, cell: usize, s: u32, window: i32) {
    // SAFETY: cell < width.
    rare_apply(unsafe { &mut *table.add(cell) }, s, window);
    if s.is_multiple_of(3) {
        // SAFETY: (cell + 1) % width < width.
        unsafe { (*table.add((cell + 1) % width)).set_raw_bits(s.wrapping_add(window as u32)) };
    }
}

// No NCBI counterpart: the forms of the branch-free loop under test.
/// Form 0 / 1: the loop before LOSAT_X_TBXHITSTEP (`x_hits_before`) without / with `PF`; 2: its
/// wasm32 form `x_hits_bf_split`; 3 / 4: `x_hits_steps` without / with `PF`; 5 / 6: `x_hits_bf`
/// with `STEP` (and `PF`) / without.
///
/// SAFETY: as `x_hits_bf` for the first `n` of `pairs`, `X_PREFETCH_AHEAD` entries after them.
#[allow(clippy::too_many_arguments)]
unsafe fn run_form<C: XCells>(
    form: usize,
    ptr: *mut DiagStruct,
    pairs: &[OffsetPair],
    n: usize,
    cells: C,
    width: usize,
    window: i32,
    wordsize: i32,
    rare: impl FnMut(usize, u32, u32),
) {
    assert!(pairs.len() >= n + X_PREFETCH_AHEAD);
    let pairs = pairs.as_ptr();
    // SAFETY: the caller's bounds, passed on.
    unsafe {
        match form {
            0 => x_hits_before::<C, false>(ptr, pairs, n, cells, width, window, wordsize, rare),
            1 => x_hits_before::<C, true>(ptr, pairs, n, cells, width, window, wordsize, rare),
            2 => x_hits_bf_split(ptr, pairs, n, cells, width, window, wordsize, rare),
            3 => x_hits_steps::<C, false>(ptr, pairs, n, cells, width, window, wordsize, rare),
            4 => x_hits_steps::<C, true>(ptr, pairs, n, cells, width, window, wordsize, rare),
            5 => x_hits_bf::<C, true, true>(ptr, pairs, n, cells, width, window, wordsize, rare),
            _ => x_hits_bf::<C, false, false>(ptr, pairs, n, cells, width, window, wordsize, rare),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-552 (module doc)
// `x_hits_bf` against the reference's statements on random tables and hits over the range of
// `fits`, with the last_hit values on both sides of each boundary of the two tests.
#[test]
fn x_hits_bf_takes_the_reference_outcomes() {
    let mut rng = Rng(0x0b1f_5eed_7777_0001);
    // too far, overlap, flag 1, candidate; s - L at 0, -1, wordsize - window, one below it.
    let mut taken = [0u64; 4];
    let mut edges = [0u64; 4];
    // Batches of no hit, batches whose last hit goes to the per-hit body, and two hits in a
    // row that go to it.
    let mut shapes = [0u64; 3];
    for round in 0..3000 {
        let wordsize = [2, 3, 3, 5][round % 4];
        let window = match round % 5 {
            0 => wordsize + 1,
            1 | 4 => 40,
            2 => 64,
            _ => 1 << 20,
        };
        // Every seventh round the longest chunk `fits` takes; every thirteenth the widest window.
        let chunk_len: usize = if round % 7 == 0 {
            ((X_HITBF_LIMIT - 1) - window as i64) as usize
        } else {
            1 + rng.below(5000) as usize
        };
        let window = if round % 13 == 6 {
            ((X_HITBF_LIMIT - 1) as usize - chunk_len) as i32
        } else {
            window
        };
        assert!(fits(chunk_len, window, wordsize));
        assert!(!fits(chunk_len + 1, window, wordsize) || round % 7 != 0);
        let width = 1 + rng.below(48) as usize;
        // The unit's lowest query offset q_min (below 2^31, as a frame's offsets are; every fifth
        // round below the chunk length, so that lo is negative) and lo = q_min - (chunk_len - 1),
        // as `scan_unit` sets it; `lo_exact` is the value `XCells::lowest_diagonal` takes.
        let q_min: i64 = if round % 5 == 3 {
            rng.below(chunk_len as u64) as i64
        } else {
            rng.below((1u64 << 31) - 64 - chunk_len as u64) as i64
        };
        let lo_exact = q_min - (chunk_len as i64 - 1);
        let lo = lo_exact as u32;
        assert_eq!(XUnmasked { lo }.lowest_diagonal() as i64, lo_exact);
        let n = rng.below(400) as usize;
        let mut table = vec![DiagStruct::default(); width];
        let mut touched = vec![false; width];
        let mut pairs: Vec<OffsetPair> = Vec::with_capacity(n + X_PREFETCH_AHEAD);
        for _ in 0..n {
            let cell = rng.below(width as u64) as usize;
            let s = match rng.below(8) {
                0 => 0,
                1 => chunk_len as u32 - 1,
                _ => rng.below(chunk_len as u64) as u32,
            };
            // A query offset is not negative: s at least -(lo + cell), which is below chunk_len.
            let s = s.max((-(lo_exact + cell as i64)).max(0) as u32);
            if !touched[cell] {
                touched[cell] = true;
                // The first hit on a cell: a start value on an edge of a test, or any value a
                // unit stores.
                let delta: i64 = match rng.below(9) {
                    0 => 0,
                    1 => -1,
                    2 => 1,
                    3 => (wordsize - window) as i64,
                    4 => (wordsize - window) as i64 - 1,
                    5 => (wordsize - window) as i64 + 1,
                    _ => i64::MAX,
                };
                let last_hit = s as i64 - delta;
                let c = &mut table[cell];
                if delta != i64::MAX && last_hit >= 0 {
                    c.set_raw_bits(last_hit as u32);
                } else if rng.below(4) == 0 {
                    c.set_raw_bits(0);
                } else {
                    let v = rng.below(chunk_len as u64 + window as u64) as u32;
                    c.set_raw_bits(v | if rng.below(5) == 0 { 0x8000_0000 } else { 0 });
                }
            }
            let q = cell as i64 + s as i64 + lo_exact;
            assert!((0..1 << 31).contains(&q), "a query offset in [0, 2^31)");
            let q = q as u32;
            assert_eq!(XUnmasked { lo }.cell(q, s), cell);
            pairs.push(OffsetPair { q_off: q, s_off: s });
        }
        // Stale entries after the batch, as the vector keeps them.
        for _ in 0..X_PREFETCH_AHEAD {
            pairs.push(OffsetPair {
                q_off: rng.next() as u32,
                s_off: rng.next() as u32,
            });
        }

        // The reference's statements (diag_offset = window).
        let mut reference = table.clone();
        let mut reference_log = Vec::new();
        let mut last_rare = None;
        shapes[0] += u64::from(n == 0);
        let mut to_rare = |i: usize| {
            shapes[1] += u64::from(i + 1 == n);
            shapes[2] += u64::from(i > 0 && last_rare == Some(i - 1));
            last_rare = Some(i);
        };
        for (i, pair) in pairs[..n].iter().enumerate() {
            let (q, s) = (pair.q_off, pair.s_off);
            let cell = q.wrapping_sub(s).wrapping_sub(lo) as usize;
            let c = &mut reference[cell];
            if c.flag() != 0 {
                taken[2] += 1;
                to_rare(i);
                reference_log.push((cell, q, s));
                // SAFETY: cell < width, the table's length.
                unsafe { rare_and_other(reference.as_mut_ptr(), width, cell, s, window) };
                continue;
            }
            let d = s as i64 - c.last_hit() as i64;
            for (k, edge) in [
                0,
                -1,
                (wordsize - window) as i64,
                (wordsize - window) as i64 - 1,
            ]
            .into_iter()
            .enumerate()
            {
                edges[k] += u64::from(d == edge);
            }
            let last_hit = c.last_hit() - window;
            let diff = s.wrapping_sub(last_hit as u32) as i32;
            if diff >= window {
                taken[0] += 1;
                c.set_last_hit(s.wrapping_add(window as u32) as i32);
                continue;
            }
            if diff < wordsize {
                taken[1] += 1;
                continue;
            }
            taken[3] += 1;
            to_rare(i);
            reference_log.push((cell, q, s));
            // SAFETY: cell < width, the table's length.
            unsafe { rare_and_other(reference.as_mut_ptr(), width, cell, s, window) };
        }

        // The loop before LOSAT_X_TBXHITSTEP with and without `PF`, its wasm32 form
        // (`x_hits_bf_split`), `x_hits_steps` with and without `PF`, and `x_hits_bf` with and
        // without `STEP`. LOSAT_X_TBXUNITMASK: the `XMasked` instances too, with the mask of all
        // ones (the same cells as `XUnmasked`).
        for form in 0..7 {
            for masked in [false, true] {
                let mut got = table.clone();
                let mut log = Vec::new();
                let ptr = got.as_mut_ptr();
                let rare = |cell: usize, q: u32, s: u32| {
                    log.push((cell, q, s));
                    // SAFETY: cell < width, the table's length (the loop's bound).
                    unsafe { rare_and_other(ptr, width, cell, s, window) };
                };
                // SAFETY: n entries and X_PREFETCH_AHEAD more; every cell is below width; `fits`
                // holds, every flag-0 value is below chunk_len + window, and every hit's q - s -
                // lo_exact in exact integers is its cell.
                unsafe {
                    if masked {
                        let all = XMasked { lo, mask: u32::MAX };
                        run_form(form, ptr, &pairs, n, all, width, window, wordsize, rare)
                    } else {
                        let cells = XUnmasked { lo };
                        run_form(form, ptr, &pairs, n, cells, width, window, wordsize, rare)
                    }
                };
                assert_eq!(
                    log, reference_log,
                    "round {round} form {form} masked {masked}: rare calls"
                );
                assert_eq!(
                    got, reference,
                    "round {round} form {form} masked {masked}: table"
                );
            }
        }
    }
    assert!(taken.iter().all(|&t| t > 1000), "outcomes {taken:?}");
    assert!(edges.iter().all(|&e| e > 100), "edges {edges:?}");
    assert!(shapes.iter().all(|&e| e > 0), "batch shapes {shapes:?}");
}

// No NCBI counterpart: compares the branch-free scans with the scans before the switch on the
// same units (module doc), and the batched ones (LOSAT_X_TBXEXTBATCH) with batches of 1, 3 and
// 16 candidates.
#[test]
fn bf_scans_equal_the_reference_scans() {
    let alphabet: Vec<u8> = "ACDEKLMNPQRSTVWY"
        .bytes()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect();
    let narrow: Vec<u8> = "AKLE"
        .bytes()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect();
    let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let mut rng = Rng(0x41b0_f00d_2026_1019);
    let mut outcomes = [0u64; 7];
    let mut batched_units = 0usize;
    for case in 0..16 {
        let abc = if case % 3 == 2 { &narrow } else { &alphabet };
        let q_len = 150 + rng.below(350) as usize;
        let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
            .iter()
            .map(|&f| frame(f, &residues(&mut rng, q_len, abc, &[])))
            .collect();
        let source: Vec<u8> = frames
            .iter()
            .flat_map(|f| f.aa_seq[1..f.aa_len + 1].to_vec())
            .collect();
        let s_len = 300 + rng.below(500) as usize;
        // Every fourth case scans the query against itself (long extensions, flag-1 cells hit
        // again inside and past their extensions).
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
        // The default window, the smallest `fits` takes, a wide one, and the smallest of the
        // batched scan.
        let window = [40, 4, 64, 3 + x_tbx_extbatch::BATCH as i32][case % 4];
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
                ncbi_length: (s_frame.aa_len - offset) as i32,
                overlap_past_list: false,
            };
            assert!(fits(chunk.length, window, 3));
            for (g, group) in pp.groups.iter().enumerate() {
                if group.lookup.longest_chain <= 0 {
                    continue;
                }
                let unit = pp.x_unit(&inp, f, 0, chunk, g);
                let mut count = XCount::default();
                let mut reference: Vec<InitHSP> = Vec::new();
                pp.scan_unit::<_, _, false>(&unit, &inp, &mut scratch, &mut reference, &mut count);
                let reference_diag = scratch.diag.clone();
                for (o, c) in outcomes.iter_mut().zip(count.0) {
                    *o += c;
                }
                // LOSAT_X_TBXHITSTEP: both forms of the loop (`STEP`).
                for (pf, step) in [(false, true), (true, true), (false, false), (true, false)] {
                    let mut got: Vec<InitHSP> = Vec::new();
                    match (pf, step) {
                        (true, true) => {
                            pp.scan_unit_bf::<_, true, true, 0>(&unit, &inp, &mut scratch, &mut got)
                        }
                        (false, true) => pp.scan_unit_bf::<_, false, true, 0>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                        ),
                        (true, false) => pp.scan_unit_bf::<_, true, false, 0>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                        ),
                        (false, false) => pp.scan_unit_bf::<_, false, false, 0>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                        ),
                    }
                    assert_eq!(
                        got, reference,
                        "case {case} frame {f} group {g} pf {pf} step {step}"
                    );
                    assert!(
                        scratch.diag == reference_diag,
                        "case {case} frame {f} group {g} pf {pf} step {step}: table"
                    );
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
                pp.scan_unit_bf::<_, false, true, 0>(&unit, &inp, &mut scratch, &mut compact);
                assert_eq!(
                    compact, compact_ref,
                    "case {case} frame {f} group {g}: compact"
                );
                // The batched scan's windows leave room for its pending marks.
                if window < 3 + x_tbx_extbatch::BATCH as i32 {
                    continue;
                }
                batched_units += 1;
                for len in [1, 3, 16] {
                    let mut got: Vec<InitHSP> = Vec::new();
                    if len == 3 {
                        pp.scan_unit_batched_hits::<_, true, true, true, 0>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            len,
                        );
                    } else if len == 16 {
                        // LOSAT_X_TBXHITSTEP=0's loop.
                        pp.scan_unit_batched_hits::<_, false, true, false, 0>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            len,
                        );
                    } else {
                        pp.scan_unit_batched_hits::<_, false, true, true, 0>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            len,
                        );
                    }
                    assert_eq!(
                        got, reference,
                        "case {case} frame {f} group {g} batch {len}"
                    );
                    assert!(
                        scratch.diag == reference_diag,
                        "case {case} frame {f} group {g} batch {len}: table"
                    );
                }
            }
        }
    }
    // Every outcome of the reference loop but the context check (one query: the unit's frames
    // are contexts of one group only where they are apart) was reached.
    let [masked, reset, too_far, overlap, _context, extended, saved] = outcomes;
    assert!(
        masked > 100 && reset > 100 && too_far > 10_000 && overlap > 1000,
        "outcomes {outcomes:?}"
    );
    assert!(extended > 1000 && saved > 100, "outcomes {outcomes:?}");
    assert!(batched_units > 20, "batched units {batched_units}");
}

// No NCBI counterpart: the rule (module doc) at its boundaries.
#[test]
fn fits_at_its_boundaries_and_units_outside_run_the_reference() {
    let top = (X_HITBF_LIMIT - 1) as usize;
    assert!(fits(top - 40, 40, 3));
    assert!(!fits(top - 39, 40, 3));
    assert!(fits(0, 40, 3));
    assert!(fits(100, 4, 3));
    assert!(!fits(100, 3, 3));
    assert!(!fits(100, 0, 3));
    assert!(!fits(100, 40, 0));
    assert!(!fits(usize::MAX, 40, 3));
    assert!(!fits(100, i32::MAX, 3));
    for on in [false, true] {
        assert_eq!(unit_bf(on, 100, 40, 3), on);
        assert!(!unit_bf(on, top - 39, 40, 3));
        assert!(!unit_bf(on, 100, 3, 3));
    }
}
