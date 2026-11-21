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

fn aa(s: &str) -> Vec<u8> {
    s.bytes()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect()
}

// No NCBI counterpart: the rule's arithmetic (module doc).
#[test]
fn the_partition_cuts_only_groups_wider_than_a_slice() {
    assert_eq!(x_partition(&[(1, 101), (102, 150)], 100, 0), None);
    let parts = x_partition(&[(1, 102), (103, 150)], 100, 0).expect("a cut");
    assert_eq!(parts, vec![vec![1, 51], vec![103]]);
    // Equal slices, the last ending at the group's end.
    let parts = x_partition(&[(10, 1010)], 300, 0).expect("a cut");
    assert_eq!(parts, vec![vec![10, 260, 510, 760]]);
    // At most X_QSLICE_MAX_SLICES slices; strictly ascending starts.
    let parts = x_partition(&[(0, 1_000_000)], 2, 0).expect("a cut");
    assert_eq!(parts[0].len() as u64, X_QSLICE_MAX_SLICES);
    assert!(parts[0].windows(2).all(|w| w[0] < w[1]));
    // A group shorter than the cap of slices: one letter each.
    let parts = x_partition(&[(5, 9)], 1, 0).expect("a cut");
    assert_eq!(parts, vec![vec![5, 6, 7, 8]]);
    // An empty group stays one slice.
    assert_eq!(x_partition(&[(7, 7)], 2, 0), None);
    // The built-in rule's bound: only the groups of at least `min_len` letters are cut.
    assert_eq!(x_partition(&[(0, 1000), (1000, 1500)], 100, 1001), None);
    let parts = x_partition(&[(0, 1000), (1000, 1500)], 100, 600).expect("a cut");
    assert_eq!(parts[0].len(), 10);
    assert_eq!(parts[1], vec![1000]);
    assert_eq!(
        X_QSLICE_MIN_SPAN_THREADS.div_ceil(8),
        768 << 10,
        "8 threads: groups of 768 Ki letters or more"
    );
}

// No NCBI counterpart: the cut of the lookup locations (module doc, item 5): every word start of
// every location goes to exactly one piece, in the slice that holds it, in list order, for
// words of 2, 3 and 4 letters (`-word_size`).
#[test]
fn slice_locations_keep_every_word_once_in_its_slice() {
    let mut rng = Rng(0x51ce_0000_1234_abcd);
    for round in 0..3000 {
        let wl = [3i32, 2, 4][round % 3];
        let a = rng.below(50) as i32;
        let b = a + 1 + rng.below(400) as i32;
        let k = 1 + rng.below(8) as usize;
        let mut parts: Vec<i32> = vec![a];
        let mut cuts: Vec<i32> = (1..k)
            .map(|_| a + 1 + rng.below((b - a) as u64) as i32)
            .collect();
        cuts.sort_unstable();
        cuts.dedup();
        cuts.retain(|&c| c > a && c < b);
        parts.extend(cuts);
        let mut locations = Vec::new();
        let mut at = a;
        while at < b {
            let from = at + rng.below(4) as i32;
            let to = (from + rng.below(30) as i32).min(b - 1);
            if from <= to {
                locations.push((from, to));
            }
            at = to + 1 + rng.below(3) as i32;
        }
        let pieces = x_slice_locations(&locations, &parts, wl);
        let starts = |list: &[(i32, i32)]| -> Vec<i32> {
            list.iter()
                .flat_map(|&(from, to)| from..=(to - wl + 1))
                .collect()
        };
        let mut all = Vec::new();
        for (j, list) in pieces.iter().enumerate() {
            let lo = parts[j];
            let hi = parts.get(j + 1).copied().unwrap_or(i32::MAX);
            for &(from, to) in list {
                assert!(to - from + 1 >= wl, "a piece holds a word");
                assert!(
                    locations.iter().any(|&(f, t)| f <= from && to <= t),
                    "a piece lies in a location"
                );
                for s in from..=(to - wl + 1) {
                    assert!(
                        j == 0 && s < lo || (lo <= s && s < hi),
                        "word {s} in slice {j}"
                    );
                }
            }
            let mut ordered = starts(list);
            let sorted = ordered.windows(2).all(|w| w[0] < w[1]);
            assert!(sorted, "slice {j}: list order is offset order");
            all.append(&mut ordered);
        }
        all.sort_unstable();
        assert_eq!(all, starts(&locations));
    }
}

// No NCBI counterpart: one instance of the scan in slices by the bits of `variant`: batched
// extensions (2), the prefetch (1), the branch-free hit loop of LOSAT_X_TBXHITBF (4); `SENT`:
// the extensions of LOSAT_X_TBXUNGAPSENT (for a unit `x_unit_sent` accepts).
fn scan_sliced_variant<S: XInitSink, M: TbxScore, const SENT: usize>(
    pp: &XPairPar<M>,
    variant: usize,
    sliced: &XSliced,
    unit: &XUnit,
    inp: &XSubject<'_>,
    scratch: &mut XScratch,
    got: &mut S,
) {
    match variant {
        0 => pp
            .scan_unit_sliced_hits::<_, false, false, false, SENT>(sliced, unit, inp, scratch, got),
        1 => {
            pp.scan_unit_sliced_hits::<_, false, true, false, SENT>(sliced, unit, inp, scratch, got)
        }
        2 => {
            pp.scan_unit_sliced_hits::<_, true, false, false, SENT>(sliced, unit, inp, scratch, got)
        }
        3 => {
            pp.scan_unit_sliced_hits::<_, true, true, false, SENT>(sliced, unit, inp, scratch, got)
        }
        4 => {
            pp.scan_unit_sliced_hits::<_, false, false, true, SENT>(sliced, unit, inp, scratch, got)
        }
        5 => {
            pp.scan_unit_sliced_hits::<_, false, true, true, SENT>(sliced, unit, inp, scratch, got)
        }
        6 => {
            pp.scan_unit_sliced_hits::<_, true, false, true, SENT>(sliced, unit, inp, scratch, got)
        }
        _ => pp.scan_unit_sliced_hits::<_, true, true, true, SENT>(sliced, unit, inp, scratch, got),
    }
}

// No NCBI counterpart: compares the scan in slices with `scan_unit` on the same units (module
// doc, items 1-4), for slices built directly and filtered (item 5), with and without the
// batched extensions, the prefetch and the branch-free hit loop (LOSAT_X_TBXHITBF).
#[test]
fn slices_equal_scan_unit() {
    slices_check(Blosum62Score, 0x5eed_51ce_0123_4567);
}

// No NCBI counterpart: `slices_equal_scan_unit` for a search of 2-letter words (`-word_size
// 2`, the other matrices' instance; SC10b PR-2): the slices of the direct builder at that word
// length, the sentinel walks of 2-letter words (`SENT` = 2).
#[test]
fn slices_of_2_letter_words_equal_scan_unit() {
    use crate::algorithm::tblastx::scoring::MatrixScore;
    for (i, matrix) in [ScoringMatrix::Blosum62, ScoringMatrix::Pam30]
        .into_iter()
        .enumerate()
    {
        let scores = MatrixScore::with_word(matrix, 2);
        slices_check(scores, 0x2b17_51ce_0000_0001 + i as u64);
    }
}

/// The checks of `slices_equal_scan_unit` for the search of `sc` (its matrix and word length).
fn slices_check<M: TbxScore>(sc: M, seed: u64) {
    let alphabet = aa("ACDEKLMNPQRSTVWY");
    let narrow = aa("AKLE");
    let (matrix, word) = (sc.matrix(), sc.word());
    let params = lookup_protein_params_ungapped(matrix);
    let serial = crate::utils::threading::SearchPool::x_serial();
    let mut rng = Rng(seed);
    let mut records = 0usize;
    let mut cut_groups = 0usize;
    // LOSAT_X_TBXUNGAPSENT: the units, those whose chunk is the whole frame, and those that took
    // the sentinel walks.
    let (mut units, mut whole_units, mut sent_units) = (0usize, 0usize, 0usize);
    let on = x_tbx_sent_contexts(std::iter::empty());
    for case in 0..12 {
        let abc = if case % 3 == 2 { &narrow } else { &alphabet };
        let q_len = 200 + rng.below(500) as usize;
        let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
            .iter()
            .map(|&f| frame(f, &residues(&mut rng, q_len, abc, &[])))
            .collect();
        let source: Vec<u8> = frames
            .iter()
            .flat_map(|f| f.aa_seq[1..f.aa_len + 1].to_vec())
            .collect();
        let s_len = 300 + rng.below(600) as usize;
        let s_frames: Vec<QueryFrame> = if case % 4 == 3 {
            frames.clone()
        } else {
            [1i8, 2, 3, -1, -2, -3]
                .iter()
                .map(|&f| frame(f, &residues(&mut rng, s_len, abc, &source)))
                .collect()
        };
        let queries = [frames.clone()];
        let threshold = 11;
        let (lookup, contexts) = crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile(
            &queries, threshold, &params, true, matrix, word,
        );
        let starts = group_starts(&contexts, X_MAX_GROUPS);
        // Every fourth case cuts no group (one slice each, as a group of a cut batch that is
        // not cut).
        let spans = x_group_spans(&contexts, &starts);
        let width = [7usize, 60, 150, 1000][case % 4];
        let parts = x_partition(&spans, width, 0)
            .unwrap_or_else(|| spans.iter().map(|&(a, _)| vec![a]).collect());
        let filtered = build_sliced(&lookup, &parts);
        let q = XPairLutQuery::x_prepare_for(&queries, threshold, true, matrix, word)
            .expect("direct builder");
        let direct = lut_direct_sliced(&q, &starts, &parts, &serial, false).expect("groups");
        lut_shadow_compare(Some(&direct), Some(&filtered));
        cut_groups += direct.iter().filter(|g| g.slices.len() > 1).count();
        let groups = build_groups(&lookup, &contexts, X_MAX_GROUPS);
        // The stub's longest_chain over the slices is the full table's.
        let (stub, _) = q.x_into_stub(tables(&direct));
        assert_eq!(stub.longest_chain, lookup.longest_chain);
        for (g, s) in groups.iter().zip(&direct) {
            assert_eq!((g.q_min, g.q_max), (s.q_min, s.q_max), "group bounds");
            assert_eq!(g.lookup.longest_chain > 0, s.has_entries());
        }
        let cutoff_scores: Vec<i32> = (0..contexts.len())
            .map(|i| 10 + (i as i32 % 4) * 3)
            .collect();
        let x_dropoff_per_context: Vec<i32> = (0..contexts.len())
            .map(|i| [7, 16, 30, 200][i % 4])
            .collect();
        let pp = XPairPar {
            groups,
            stub: false,
            // As `prepare` decides it.
            sent: sc.x_sent_walks()
                && x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
            sc,
            sliced: Some(direct),
            mask_rule: None,
            extents: Vec::new(),
            ctx_of: x_ctx_table(&lookup, true),
        };
        let window = [40, 3 + BATCH as i32, 64][case % 3];
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
        let mut scratch = XScratch::new(false);
        for (f, s_frame) in s_frames.iter().enumerate() {
            let offset = if case % 2 == 1 { s_frame.aa_len / 3 } else { 0 };
            let chunk = SubjectChunk {
                offset,
                length: s_frame.aa_len - offset,
                overlap: 0,
            };
            for g in 0..pp.groups.len() {
                if pp.groups[g].lookup.longest_chain <= 0 {
                    continue;
                }
                let unit = XUnit { f, k: 0, chunk, g };
                let mut reference: Vec<InitHSP> = Vec::new();
                pp.scan_unit::<_, _, false>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut reference,
                    &mut XNoProbe,
                );
                let reference_diag = scratch.diag.clone();
                reference.x_sort();
                for sliced in [&pp.sliced.as_ref().unwrap()[g], &filtered[g]] {
                    let mut plain: Vec<Vec<InitHSP>> = Vec::new();
                    // 8..16: the same with LOSAT_X_TBXUNGAPSENT, on the units that take it.
                    let sent = pp.x_unit_sent(&unit, &inp);
                    units += 1;
                    whole_units += usize::from(chunk.offset == 0);
                    sent_units += usize::from(sent);
                    for variant in 0..16 {
                        if variant >= 8 && !sent {
                            continue;
                        }
                        let mut got: Vec<InitHSP> = Vec::new();
                        if variant < 8 {
                            scan_sliced_variant::<_, _, 0>(
                                &pp,
                                variant,
                                sliced,
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                            );
                        } else if word == 2 {
                            scan_sliced_variant::<_, _, 2>(
                                &pp,
                                variant & 7,
                                sliced,
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                            );
                        } else {
                            scan_sliced_variant::<_, _, 3>(
                                &pp,
                                variant & 7,
                                sliced,
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                            );
                        }
                        assert!(
                            scratch.diag == reference_diag,
                            "case {case} frame {f} group {g} variant {variant}: table"
                        );
                        if variant < 4 {
                            plain.push(got.clone());
                        } else {
                            // LOSAT_X_TBXHITBF: the same records in the same order.
                            assert_eq!(
                                got,
                                plain[variant & 3],
                                "case {case} frame {f} group {g} variant {variant}: order"
                            );
                        }
                        got.x_sort();
                        assert_eq!(
                            got, reference,
                            "case {case} frame {f} group {g} variant {variant}: records"
                        );
                    }
                }
                records += reference.len();
                // The compact records, sorted, are today's too.
                let mut compact_ref: Vec<XUnitHsp> = Vec::new();
                pp.scan_unit::<_, _, false>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut compact_ref,
                    &mut XNoProbe,
                );
                let mut compact: Vec<XUnitHsp> = Vec::new();
                pp.scan_unit_sliced::<_, true, false>(
                    &pp.sliced.as_ref().unwrap()[g],
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut compact,
                );
                let mut compact_bf: Vec<XUnitHsp> = Vec::new();
                pp.scan_unit_sliced_hits::<_, true, false, true, 0>(
                    &pp.sliced.as_ref().unwrap()[g],
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut compact_bf,
                );
                assert_eq!(
                    compact_bf, compact,
                    "case {case} frame {f} group {g}: compact, branch-free"
                );
                compact.x_sort();
                compact_ref.x_sort();
                assert_eq!(
                    compact, compact_ref,
                    "case {case} frame {f} group {g}: compact"
                );
            }
        }
    }
    assert!(cut_groups > 40, "cut groups {cut_groups}");
    assert!(records > 500, "records {records}");
    // The units of a whole frame take the walks when the switch is on, the others (the odd
    // cases' chunks start inside the frame) never.
    assert!(0 < whole_units && whole_units < units);
    assert_eq!(sent_units, if on { whole_units } else { 0 });
}
