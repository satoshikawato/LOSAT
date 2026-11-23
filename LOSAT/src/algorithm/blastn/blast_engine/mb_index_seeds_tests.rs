use super::*;

#[test]
fn a_chunks_number_counts_ncbis_chunk_offsets() {
    let step = MAX_DBSEQ_LEN - 100;
    assert_eq!(chunk_number(0, 100), 0);
    assert_eq!(chunk_number(step, 100), 1);
    assert_eq!(chunk_number(2 * step, 100), 2);
    // An overlap whose step is not a whole byte: each offset rounds down to one.
    let step = MAX_DBSEQ_LEN - 101;
    assert_eq!(chunk_number(step - step % 4, 101), 1);
    let next = step - step % 4 + step;
    assert_eq!(chunk_number(next - next % 4, 101), 2);
}

#[test]
fn the_diagonal_hash_keeps_ncbis_entries() {
    let mut hash = IrHash::new();
    let d1 = 0x1000_0000u32 + 5;
    let d2 = d1 + IR_HASH_SIZE;
    // An empty key gives its first entry.
    let e = hash.locate(d1);
    assert_eq!((e.diag, e.qend), (0, 0));
    e.diag = d1;
    e.qend = 10;
    // Another diagonal of the key: a new entry after the first.
    let e = hash.locate(d2);
    assert_eq!((e.diag, e.qend), (d2, 0));
    e.qend = 20;
    assert_eq!(hash.pool.len(), 1);
    // Found after the first: its values move to the first entry.
    let e = hash.locate(d2);
    assert_eq!((e.diag, e.qend), (d2, 20));
    let key = (d1 % IR_HASH_SIZE) as usize;
    assert_eq!((hash.heads[key].diag, hash.heads[key].qend), (d2, 20));
    assert_eq!((hash.pool[0].diag, hash.pool[0].qend), (d1, 10));
    let e = hash.locate(d1);
    assert_eq!((e.diag, e.qend), (d1, 10));
    assert_eq!(hash.pool.len(), 1);
    // A first entry whose end is 0 is taken by any diagonal of its key.
    let d3 = 0x1000_0000u32 + 6;
    let e = hash.locate(d3);
    assert_eq!((e.diag, e.qend), (0, 0));
    hash.clear();
    assert_eq!(hash.heads[key].next, 0);
    assert_eq!((hash.heads[key].diag, hash.heads[key].qend), (0, 0));
    assert!(hash.pool.is_empty() && hash.used.is_empty());
}

#[test]
fn the_lookup_segments_are_the_unmasked_intervals_of_valid_contexts() {
    let context = |query_offset, len, masks: &[(usize, usize)]| QueryContext {
        query_idx: 0,
        frame: 1,
        query_offset,
        len,
        masks: masks
            .iter()
            .map(|&(start, end)| crate::utils::dust::MaskedInterval { start, end })
            .collect(),
    };
    let contexts = [
        context(0, 50, &[(10, 20)]),
        context(51, 50, &[]),
        context(102, 30, &[(0, 30)]),
        context(133, 20, &[(15, 20)]),
    ];
    assert_eq!(
        lookup_segments(&contexts, &[true, false, true, true]),
        [(0, 10), (20, 50), (133, 148)]
    );
}

/// A volume of `n_oids` sequences from `start_oid`, with or without an index.
fn volume(start_oid: u32, n_oids: u32, has_index: bool) -> IndexedVolume {
    IndexedVolume {
        path: String::new(),
        start_oid,
        n_oids,
        has_index,
    }
}

/// The plan as letters: S searched (K: and kept), E earlier results, X scanned.
fn letters(plan: &[VolumePlan]) -> String {
    plan.iter()
        .map(|plan| match plan {
            VolumePlan::Searched { kept: false } => 'S',
            VolumePlan::Searched { kept: true } => 'K',
            VolumePlan::Earlier(_) => 'E',
            VolumePlan::Scanned => 'X',
        })
        .collect()
}

/// The counts of `held`.
fn counts(held: &[HeldResults]) -> Vec<i64> {
    held.iter().map(|held| held.ref_count).collect()
}

/// A search stores the results of the volumes it keeps (`IndexSeeds::check`).
fn keep(plan: &[VolumePlan], held: &mut [HeldResults]) {
    for (at, plan) in plan.iter().enumerate() {
        if let VolumePlan::Searched { kept: true } = plan {
            held[at].results = Some(Arc::new(VolumeSeeds::empty_for_tests(2)));
        }
    }
}

#[test]
fn the_volume_before_one_without_an_index_keeps_its_results() {
    // An index volume, a database volume without an index, an index volume.
    let volumes = [volume(0, 2, true), volume(2, 2, false), volume(4, 2, true)];
    let mut held: Vec<HeldResults> = volumes.iter().map(|_| HeldResults::default()).collect();
    let plan = plan_search(&volumes, &[true; 3], &mut held);
    assert_eq!(letters(&plan), "KXS");
    assert_eq!(counts(&held), [1, -1, 0]);
    keep(&plan, &mut held);
    // Every later search takes the first volume's earlier results.
    for later in [2, 3] {
        let plan = plan_search(&volumes, &[true; 3], &mut held);
        assert_eq!(letters(&plan), "EXS");
        assert_eq!(counts(&held), [1, -later, 0]);
    }
    // The last volume without an index: the volume before it keeps its results.
    let volumes = [volume(0, 2, true), volume(2, 2, true), volume(4, 2, false)];
    let mut held: Vec<HeldResults> = volumes.iter().map(|_| HeldResults::default()).collect();
    let plan = plan_search(&volumes, &[true; 3], &mut held);
    assert_eq!(letters(&plan), "SKX");
    keep(&plan, &mut held);
    assert_eq!(
        letters(&plan_search(&volumes, &[true; 3], &mut held)),
        "SEX"
    );
}

#[test]
fn a_search_that_reaches_no_volume_without_an_index_keeps_no_results() {
    let volumes = [volume(0, 2, true), volume(2, 2, true)];
    let mut held: Vec<HeldResults> = volumes.iter().map(|_| HeldResults::default()).collect();
    for _ in 0..2 {
        assert_eq!(letters(&plan_search(&volumes, &[true; 2], &mut held)), "SS");
        assert_eq!(counts(&held), [0, 0]);
    }
    // A database volume without an index that the search does not reach (no sequence of it in
    // the OID list) counts as passed.
    let volumes = [volume(0, 2, true), volume(2, 2, false), volume(4, 2, true)];
    let mut held: Vec<HeldResults> = volumes.iter().map(|_| HeldResults::default()).collect();
    let plan = plan_search(&volumes, &[true, false, true], &mut held);
    assert_eq!(letters(&plan), "SXS");
    assert_eq!(counts(&held), [0, -1, 0]);
}

#[test]
fn the_scan_saves_into_a_freed_list_at_its_second_initial_hit_after_earlier_results() {
    let searched = VolumePlan::Searched { kept: false };
    let earlier = VolumePlan::Earlier(Arc::new(VolumeSeeds::empty_for_tests(2)));
    let scanned = VolumePlan::Scanned;
    let freed =
        |volumes: &[(&VolumePlan, bool, u32)]| saves_into_a_freed_list(volumes.iter().copied());
    // A list moved in from a volume searched here has memory.
    assert!(!freed(&[(&searched, true, 0), (&scanned, false, 2)]));
    // A list of earlier results has none: one hit is saved, the second frees the list.
    assert!(!freed(&[(&earlier, true, 0), (&scanned, false, 1)]));
    assert!(freed(&[(&earlier, true, 0), (&scanned, false, 2)]));
    assert!(freed(&[
        (&earlier, true, 0),
        (&scanned, false, 1),
        (&scanned, false, 1)
    ]));
    assert!(!freed(&[
        (&earlier, true, 0),
        (&scanned, false, 0),
        (&scanned, false, 1)
    ]));
    // A list moved in again from a searched volume, or none taken of the earlier results.
    assert!(!freed(&[
        (&earlier, true, 0),
        (&scanned, false, 1),
        (&searched, true, 0),
        (&scanned, false, 2)
    ]));
    assert!(!freed(&[(&earlier, false, 0), (&scanned, false, 2)]));
}

/// The seeds that `seeds` (left end in the sequence, query offset) give the hard chunks of a
/// sequence of `length` letters with the hard ranges `ranges`: each chunk's offset, residual,
/// length and seeds (query offset, offset in the chunk).
fn placed(
    length: usize,
    ranges: &[(i32, i32)],
    seeds: &[(usize, u32)],
) -> Vec<(usize, usize, usize, Vec<(u32, u32)>)> {
    let mut hard = HardSeeds::chunks_of(length, ranges, 100, 0);
    for &(start, q_off) in seeds {
        hard.place(Seed { q_off, s_off: 0 }, start, ranges);
    }
    hard.chunks
        .iter()
        .map(|chunk| {
            (
                chunk.offset,
                chunk.residual,
                chunk.length,
                chunk
                    .seeds
                    .iter()
                    .map(|seed| (seed.q_off, seed.s_off))
                    .collect(),
            )
        })
        .collect()
}

/// Each seed goes to the range that holds its left end, in the order given; a range after a mask
/// starts a chunk at a whole byte before it; a seed at a masked letter is not placed.
#[test]
fn a_seed_goes_to_the_range_that_holds_it() {
    // Masks [100, 150) and [210, 230).
    let ranges = [(0, 100), (150, 210), (230, 1000)];
    let got = placed(
        1000,
        &ranges,
        &[(60, 5), (180, 9), (110, 0), (150, 3), (999, 1), (40, 2)],
    );
    assert_eq!(
        got,
        vec![
            (0, 0, 100, vec![(5, 60), (2, 40)]),
            (148, 2, 62, vec![(9, 32), (3, 2)]),
            (228, 2, 772, vec![(1, 771)]),
        ]
    );
}

/// Each seed of a range over `MAX_DBSEQ_LEN` letters goes to the last chunk that starts at or
/// before it, its overlap with the chunk before included, as the index places it.
#[test]
fn a_long_range_places_a_seed_in_the_later_chunk() {
    let step = MAX_DBSEQ_LEN - 100;
    let length = MAX_DBSEQ_LEN + 1000;
    let ranges = [(0, 40), (64, length as i32)];
    let got = placed(length, &ranges, &[(64 + step - 10, 0), (64 + step + 10, 1)]);
    assert_eq!(
        got,
        vec![
            (0, 0, 40, vec![]),
            (64, 0, MAX_DBSEQ_LEN, vec![(0, (step - 10) as u32)]),
            (64 + step, 0, length - 64 - step, vec![(1, 10)]),
        ]
    );
}

/// The chunk loop's chunk finds its seeds by offset, residual and length; chunks of hard ranges of
/// no letters at one place hold none.
#[test]
fn a_hard_chunk_finds_its_seeds() {
    // Masks [0, 9), [9, 9) and [9, 20): two empty ranges at 9.
    let ranges = [(0, 0), (9, 9), (9, 9), (20, 300)];
    let mut hard = HardSeeds::chunks_of(300, &ranges, 100, 3);
    hard.word_size = 28;
    hard.place(Seed { q_off: 7, s_off: 0 }, 30, &ranges);
    let offsets: Vec<_> = hard
        .chunks
        .iter()
        .map(|chunk| (chunk.offset, chunk.residual, chunk.length))
        .collect();
    assert_eq!(offsets, vec![(0, 0, 0), (8, 1, 1), (8, 1, 1), (20, 0, 280)]);
    assert!(hard.chunk(8, 1, 1).is_none());
    let (seeds, word_size, at) = hard.chunk(20, 0, 280).unwrap();
    assert_eq!((word_size, at), (28, 3));
    assert_eq!(
        seeds,
        &[Seed {
            q_off: 7,
            s_off: 10
        }][..]
    );
}
