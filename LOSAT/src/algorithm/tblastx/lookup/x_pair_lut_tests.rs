use super::super::tests::{frame, residues};
use super::*;

/// One context per group: the first context of each group.
fn one_group_per_context(q: &XPairLutQuery) -> Vec<usize> {
    (0..q.contexts().len()).collect()
}

/// The table of group g as the filter of the full table: every chain filtered to the query
/// offsets of the group (`bounds[g] <= q < bounds[g + 1]`), order kept, in the layout of
/// `BlastAaLookupFinalize` (the rule `x_pair_par.rs::build_groups` applies).
fn filtered(full: &BlastAaLookupTable, bounds: &[i32], g: usize) -> XPairLutGroup {
    let lo = bounds[g];
    let hi = bounds.get(g + 1).copied().unwrap_or(i32::MAX);
    let mut table = BlastAaLookupTable {
        backbone: vec![BackboneCell::default(); full.backbone.len()],
        overflow: Vec::new(),
        pv: vec![0u32; full.pv.len()],
        frame_bases: full.frame_bases.clone(),
        num_contexts: full.num_contexts,
        query_length: full.query_length,
        word_length: full.word_length,
        alphabet_size: full.alphabet_size,
        charsize: full.charsize,
        mask: full.mask,
        longest_chain: 0,
        threshold: full.threshold,
        row_max: full.row_max.clone(),
    };
    let mut q_min = u32::MAX;
    let mut q_max = 0u32;
    for idx in 0..full.backbone.len() {
        let chain: Vec<i32> = full
            .get_hits(idx)
            .iter()
            .copied()
            .filter(|&q| lo <= q && q < hi)
            .collect();
        if chain.is_empty() {
            continue;
        }
        pv_set(&mut table.pv, idx);
        table.backbone[idx].num_used = chain.len() as i32;
        if chain.len() <= AA_HITS_PER_CELL {
            table.backbone[idx].entries[..chain.len()].copy_from_slice(&chain);
        } else {
            table.backbone[idx].set_cursor(table.overflow.len());
            table.overflow.extend_from_slice(&chain);
        }
        table.longest_chain = table.longest_chain.max(chain.len() as i32);
        for &q in &chain {
            q_min = q_min.min(q as u32);
            q_max = q_max.max(q as u32);
        }
    }
    XPairLutGroup {
        q_min: if table.longest_chain > 0 { q_min } else { 0 },
        q_max,
        table,
        rank: Vec::new(),
    }
}

/// Builds the direct groups of `queries` for the group starts `starts` (None: one group per
/// context) and compares each with the filter of the reference full table, and the stub with
/// the full table. Returns (entries of the full table, groups with entries, entries in the
/// overflow arrays of the groups).
fn check(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    starts: Option<&[usize]>,
) -> (usize, usize, usize) {
    let params = crate::stats::lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (full, full_contexts) = build_ncbi_lookup(queries, threshold, &params, true);
    let q = XPairLutQuery::x_prepare(queries, threshold, true).expect("direct builder applies");
    x_assert_same_contexts(q.contexts(), &full_contexts);
    let starts: Vec<usize> = match starts {
        Some(starts) => starts.to_vec(),
        None => one_group_per_context(&q),
    };
    let bounds: Vec<i32> = starts.iter().map(|&c| full.frame_bases[c]).collect();
    let locations = q
        .x_group_locations(&starts)
        .expect("every location inside one group");
    assert_eq!(locations.len(), starts.len());
    let groups: Vec<XPairLutGroup> = locations.iter().map(|l| q.x_build_group(l, true)).collect();
    let mut with_entries = 0usize;
    let mut in_overflow = 0usize;
    for (g, direct) in groups.iter().enumerate() {
        let mut reference = filtered(&full, &bounds, g);
        // A direct table laid out compacted equals the filtered one compacted.
        if !direct.rank.is_empty() {
            let rank = crate::algorithm::tblastx::x_scan_bf::x_compact(&mut reference.table);
            assert_eq!(direct.rank, rank, "group {g}: rank differs");
        } else {
            assert!(!crate::algorithm::tblastx::x_scan_bf::x_scan_bf_choice(
                &reference.table
            ));
        }
        x_assert_same_group(
            g,
            &direct.table,
            (direct.q_min, direct.q_max),
            &reference.table,
            (reference.q_min, reference.q_max),
        );
        with_entries += usize::from(direct.table.longest_chain > 0);
        in_overflow += direct.table.overflow.len();
    }
    let (stub, contexts) = q.x_into_stub(groups.iter().map(|g| &g.table));
    x_assert_stub_scalars(&stub, &full);
    x_assert_same_contexts(&contexts, &full_contexts);
    let entries = full.backbone.iter().map(|c| c.num_used as usize).sum();
    (entries, with_entries, in_overflow)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:92-94,105-106
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:504-509,518-519
#[test]
fn direct_groups_equal_the_filtered_full_table() {
    let mut seed = 0x5eed_0001u64;
    let mut low_complexity = vec![1u8; 200]; // poly-A: cells far above AA_HITS_PER_CELL
    low_complexity.extend(residues(&mut seed, 60));
    // One query: six contexts with masks, a frame shorter than a word, an empty frame.
    let one_query = vec![vec![
        frame(
            residues(&mut seed, 500),
            vec![(10, 30), (400, 401), (498, 500)],
            1,
        ),
        frame(residues(&mut seed, 2), vec![], 2),
        frame(low_complexity, vec![(150, 199)], 3),
        frame(residues(&mut seed, 450), vec![(0, 5)], -1),
        frame(Vec::new(), vec![], -2),
        frame(residues(&mut seed, 480), vec![], -3),
    ]];
    for threshold in [0, 11, 13] {
        let (entries, with_entries, in_overflow) = check(&one_query, threshold, None);
        assert!(entries > 0);
        assert_eq!(with_entries, 4, "threshold {threshold}");
        assert!(in_overflow > 0, "the fixture must reach the overflow array");
    }
    // Every number of contexts from 1 to 6.
    for n in 1..=6usize {
        let frames: Vec<QueryFrame> = (0..n)
            .map(|c| {
                let len = 40 + 37 * c;
                frame(
                    residues(&mut seed, len),
                    if c % 2 == 0 { vec![(5, 12)] } else { vec![] },
                    [1i8, 2, 3, -1, -2, -3][c],
                )
            })
            .collect();
        for threshold in [0, 11, 13] {
            let (_, with_entries, _) = check(std::slice::from_ref(&frames), threshold, None);
            assert_eq!(with_entries, n);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-590
// A threshold of INT_MIN (NCBI's conversion of a threshold beyond Int4) makes every word a
// neighbour of every query word.
#[test]
fn direct_groups_with_every_word_a_neighbour() {
    let mut seed = 0x5eed_0002u64;
    let query = vec![vec![
        frame(residues(&mut seed, 9), vec![], 1),
        frame(residues(&mut seed, 7), vec![(2, 4)], 2),
        frame(residues(&mut seed, 2), vec![], 3),
        frame(residues(&mut seed, 8), vec![], -1),
    ]];
    let (entries, with_entries, in_overflow) = check(&query, i32::MIN, None);
    assert_eq!(with_entries, 3);
    // 7 + 1 + 6 word starts (the mask (2, 4) of the 7-letter frame leaves the locations [0, 1],
    // too short for a word, and [4, 6], one word), each in every one of the 28^3 cells.
    assert_eq!(entries, 14 * 28 * 28 * 28);
    // The cells of the groups with 7 and 6 words are in their overflow arrays.
    assert_eq!(in_overflow, 13 * 28 * 28 * 28);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1047
// Several queries: more than six contexts, taken as groups of several contexts, and a context
// that is masked from end to end.
#[test]
fn direct_groups_of_several_contexts_and_an_all_masked_context() {
    let mut seed = 0x5eed_0003u64;
    let queries: Vec<Vec<QueryFrame>> = (0..3)
        .map(|q| {
            (0..6)
                .map(|c| {
                    let len = 30 + 11 * ((q * 6 + c) % 7);
                    let masks = if q == 1 && c == 2 {
                        vec![(0, len)] // all masked
                    } else if c == 4 {
                        vec![(3, 9), (20, 22)]
                    } else {
                        vec![]
                    };
                    frame(residues(&mut seed, len), masks, [1i8, 2, 3, -1, -2, -3][c])
                })
                .collect()
        })
        .collect();
    for threshold in [0, 11, 13] {
        // 18 groups of one context: the all-masked context gives an empty group.
        let (_, with_entries, _) = check(&queries, threshold, None);
        assert_eq!(with_entries, 17);
        // Six groups of three contexts, and uneven groups.
        let (_, with_entries, _) = check(&queries, threshold, Some(&[0, 3, 6, 9, 12, 15]));
        assert_eq!(with_entries, 6);
        check(&queries, threshold, Some(&[0, 1, 7, 8, 9, 17]));
        // One group: the direct table is the full table.
        let (_, with_entries, _) = check(&queries, threshold, Some(&[0]));
        assert_eq!(with_entries, 1);
    }
    // A query whose every context is masked: no group has an entry.
    let masked: Vec<Vec<QueryFrame>> = vec![(0..6)
        .map(|c| {
            frame(
                residues(&mut seed, 40),
                vec![(0, 40)],
                [1i8, 2, 3, -1, -2, -3][c],
            )
        })
        .collect()];
    assert_eq!(check(&masked, 13, None), (0, 0, 0));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1095-1099
// A mask interval that starts past the end of its context (SEG never makes one) gives a
// location that runs into the next context: such a batch is left to the reference builder.
#[test]
fn a_location_outside_one_group_declines() {
    let mut seed = 0x5eed_0004u64;
    let queries = vec![vec![
        frame(residues(&mut seed, 20), vec![(25, 28)], 1),
        frame(residues(&mut seed, 20), vec![], 2),
    ]];
    let q = XPairLutQuery::x_prepare(&queries, 13, true).expect("direct builder applies");
    assert!(q.x_group_locations(&[0, 1]).is_none());
    // Inside one group of both contexts the location is fine, and the tables agree.
    assert!(q.x_group_locations(&[0]).is_some());
    check(&queries, 13, Some(&[0]));
}
