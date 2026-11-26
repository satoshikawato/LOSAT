use super::*;

#[test]
fn masks_before_a_part_of_a_ranged_query() {
    let part = |from, to| ChunkQuery { query: 0, from, to };
    // The query's interval starts at record letter 100: a mask of the searched letters
    // 296..298 ends (with the extra residue) 2 letters before the part 300..600's interval.
    let masks = [MaskedInterval::new(296, 298)];
    assert_eq!(masked_before_part(&masks, 0, &part(300, 600), 100), Some(2));
    // The first part keeps the lower case before the interval: record letter 97 (bit 2)
    // ends 2 letters before the interval with its extra residue; letter 99 (bit 0) reaches
    // the part's first letter, which `restrict_masks_with_letter_before` masks.
    assert_eq!(masked_before_part(&[], 0b100, &part(0, 600), 100), Some(2));
    assert_eq!(masked_before_part(&[], 0b1, &part(0, 600), 100), None);
    // Without a range nothing lies before a part's interval.
    assert_eq!(masked_before_part(&masks, 0, &part(300, 600), 0), None);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:123-138
#[test]
fn the_chunk_count_and_size_are_ncbis() {
    assert_eq!(
        calculate_num_chunks(SplitSizes::default_for(false), 1_999_799),
        (1, 1_999_799)
    );
    assert_eq!(
        calculate_num_chunks(SplitSizes::default_for(false), 1_999_800),
        (2, 999_951)
    );
    // EDL933 with -task blastn: five chunks of 1,105,770 (E2f §D).
    assert_eq!(
        calculate_num_chunks(SplitSizes::default_for(false), 5_528_445),
        (5, 1_105_770)
    );
    assert_eq!(
        calculate_num_chunks(SplitSizes::default_for(true), 5_528_445),
        (1, 5_528_445)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:612-662
#[test]
fn a_query_over_two_chunks_is_mapped_onto_its_contexts() {
    let chunks = split_query_batch(
        &[1_999_800],
        SplitSizes::default_for(false),
        QueryStrand::Both,
    )
    .expect("split");
    assert_eq!(chunks.len(), 2);
    assert_eq!(
        chunks[0].queries,
        vec![ChunkQuery {
            query: 0,
            from: 0,
            to: 999_951
        }]
    );
    assert_eq!(
        chunks[1].queries,
        vec![ChunkQuery {
            query: 0,
            from: 999_851,
            to: 1_999_800
        }]
    );
    assert_eq!(chunks[0].contexts, vec![Some(0), Some(1)]);
    assert_eq!(chunks[0].context_offsets, vec![0, 1_999_800 - 999_951]);
    assert_eq!(chunks[1].context_offsets, vec![999_851, 0]);
}

#[test]
fn parts_that_end_or_start_in_an_overlap_keep_their_starts() {
    // Chunks of 1,000,059 overlapping by 100: the first query ends 58 residues into the
    // overlap of the second and third chunks, and the second query starts there.
    let chunks = split_query_batch(
        &[1_999_976, 1_000_000],
        SplitSizes::default_for(false),
        QueryStrand::Both,
    )
    .expect("split");
    assert_eq!(chunks.len(), 3);
    assert_eq!(
        chunks[1].queries,
        vec![
            ChunkQuery {
                query: 0,
                from: 999_959,
                to: 1_999_976
            },
            ChunkQuery {
                query: 1,
                from: 0,
                to: 42
            },
        ]
    );
    assert_eq!(
        chunks[2].queries,
        vec![
            ChunkQuery {
                query: 0,
                from: 1_999_918,
                to: 1_999_976
            },
            ChunkQuery {
                query: 1,
                from: 0,
                to: 1_000_000
            },
        ]
    );
    assert_eq!(chunks[1].contexts, vec![Some(0), Some(1), Some(2), Some(3)]);
    assert_eq!(
        chunks[1].context_offsets,
        vec![999_959, 0, 0, 1_000_000 - 42]
    );
    assert_eq!(chunks[2].context_offsets, vec![1_999_918, 0, 0, 0]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:256-281
// With `-query_loc 1000001-3200000` a mask survives in a chunk only from the chunk's
// start plus the interval's start on (oracle: `-task blastn -lcase_masking` loses the
// lower-case runs at record 1.2 Mb and 2.3 Mb, keeps those at 1.05 Mb and 2.12 Mb).
#[test]
fn ranged_chunk_masks_are_restricted_in_the_searched_coordinates() {
    let offset = 1_000_000;
    let masks = [
        MaskedInterval::new(50_000, 55_000),
        MaskedInterval::new(200_000, 205_000),
        MaskedInterval::new(1_120_000, 1_125_000),
        MaskedInterval::new(1_300_000, 1_305_000),
    ];
    let first = ChunkQuery {
        query: 0,
        from: 0,
        to: 1_100_050,
    };
    assert_eq!(
        restrict_masks(&masks, &first, offset),
        vec![MaskedInterval::new(50_000, 55_001)]
    );
    let second = ChunkQuery {
        query: 0,
        from: 1_099_950,
        to: 2_200_000,
    };
    assert_eq!(
        restrict_masks(&masks, &second, offset),
        vec![MaskedInterval::new(20_050, 25_051)]
    );
    // Without a range the same masks are restricted to each part as they are.
    assert_eq!(
        restrict_masks(&masks, &second, 0),
        vec![
            MaskedInterval::new(20_050, 25_051),
            MaskedInterval::new(200_050, 205_051)
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seqlocinfo.cpp:102-116
#[test]
fn restricted_masks_grow_by_one_and_drop_one_at_the_last_residue() {
    let part = ChunkQuery {
        query: 0,
        from: 100,
        to: 200,
    };
    let masks = [
        MaskedInterval::new(0, 50),
        MaskedInterval::new(90, 110),
        MaskedInterval::new(150, 160),
        MaskedInterval::new(190, 300),
        MaskedInterval::new(199, 250),
    ];
    assert_eq!(
        restrict_masks(&masks, &part, 0),
        vec![
            MaskedInterval::new(0, 11),
            MaskedInterval::new(50, 61),
            MaskedInterval::new(90, 100),
        ]
    );
}

// Masks in the record's coordinates are restricted as `restrict_masks` restricts the same
// masks in the searched letters' coordinates; a run that ends on the letter before the
// interval (record 9 with offset 10) reaches the first part's first letter, only while that
// part's last letter (`to - 1`, which bounds the run in record coordinates) is not before the
// offset.
#[test]
fn record_masks_restrict_as_searched_masks_and_reach_the_first_letter() {
    let offset = 10;
    let searched = [
        MaskedInterval::new(0, 3),
        MaskedInterval::new(95, 120),
        MaskedInterval::new(185, 190),
    ];
    let record: Vec<_> = searched
        .iter()
        .map(|mask| MaskedInterval::new(mask.start + offset, mask.end + offset))
        .collect();
    for (from, to) in [(0, 100), (90, 200), (180, 200), (0, 11), (0, 12)] {
        let part = ChunkQuery { query: 0, from, to };
        assert_eq!(
            restrict_record_masks(&record, &part, offset),
            restrict_masks(&searched, &part, offset),
            "{from}..{to}"
        );
    }
    let before = [MaskedInterval::new(4, 10)];
    let part = |to| ChunkQuery {
        query: 0,
        from: 0,
        to,
    };
    assert_eq!(
        restrict_record_masks(&before, &part(100), offset),
        vec![MaskedInterval::new(0, 1)]
    );
    assert_eq!(
        restrict_record_masks(&before, &part(11), offset),
        vec![MaskedInterval::new(0, 1)]
    );
    assert!(restrict_record_masks(&before, &part(10), offset).is_empty());
    let later = ChunkQuery {
        query: 0,
        from: 90,
        to: 200,
    };
    assert!(restrict_record_masks(&before, &later, offset).is_empty());
    let two_before = [MaskedInterval::new(4, 9)];
    assert!(restrict_record_masks(&two_before, &part(100), offset).is_empty());
}

// The masks of a part from the runs of the searched letters, with the letter before the
// interval when it is lower case, mask the residues that the runs of the whole record mask
// (`restrict_record_masks`, NCBI's), for every lower-case pattern of a 10-letter record, four
// intervals and every part; the searched letters' runs alone miss some.
#[test]
fn the_letter_before_the_interval_gives_the_whole_records_masks() {
    use super::super::coordination::collect_lowercase_masks;
    let residues = |masks: &[MaskedInterval]| -> std::collections::BTreeSet<usize> {
        masks.iter().flat_map(|mask| mask.start..mask.end).collect()
    };
    let mut missed = 0;
    for pattern in 0u32..1 << 10 {
        let record: Vec<u8> = (0..10)
            .map(|letter| {
                if pattern >> letter & 1 == 1 {
                    b'a'
                } else {
                    b'A'
                }
            })
            .collect();
        let runs = collect_lowercase_masks(&record);
        for (offset, end) in [(1, 10), (2, 8), (3, 10), (4, 9)] {
            let searched = collect_lowercase_masks(&record[offset..end]);
            let before = record[offset - 1].is_ascii_lowercase();
            for from in 0..end - offset {
                for to in from + 1..=end - offset {
                    let part = ChunkQuery { query: 0, from, to };
                    let whole = residues(&restrict_record_masks(&runs, &part, offset));
                    let masks = if before {
                        restrict_masks_with_letter_before(&searched, &part, offset)
                    } else {
                        restrict_masks(&searched, &part, offset)
                    };
                    assert_eq!(
                        residues(&masks),
                        whole,
                        "{} {offset}..{end} {from}..{to}",
                        String::from_utf8_lossy(&record)
                    );
                    missed +=
                        usize::from(residues(&restrict_masks(&searched, &part, offset)) != whole);
                }
            }
        }
    }
    assert!(missed > 0);
}

// CHUNK_SIZE and OVERLAP_CHUNK_SIZE as NCBI reads them (split_query_aux_priv.cpp:100-138).
#[test]
fn environment_sizes_follow_ncbis_size_t_arithmetic() {
    let sizes = SplitSizes::new(false, Some(500), None);
    assert_eq!(calculate_num_chunks(sizes, 100_000), (250, 500));
    let sizes = SplitSizes::new(false, Some(40_000), Some(0));
    assert_eq!(calculate_num_chunks(sizes, 100_000), (2, 50_001));
    // An overlap at least the chunk size, a negative overlap (2^64 - 1) and a negative
    // chunk size (2^64 - 5) never split.
    assert_eq!(
        calculate_num_chunks(SplitSizes::new(false, Some(150), Some(150)), 100_000),
        (1, 100_000)
    );
    assert_eq!(
        calculate_num_chunks(SplitSizes::new(false, Some(40_000), Some(-1)), 100_000),
        (1, 100_000)
    );
    assert_eq!(
        calculate_num_chunks(SplitSizes::new(false, Some(-5), None), 100_000),
        (1, 100_000)
    );
    // The CBatchSizeMixer maximum: size_t chunk size - 1000 as Int4.
    assert_eq!(
        SplitSizes::default_for(false).mixer_max_batch_size(),
        999_000
    );
    assert_eq!(
        SplitSizes::default_for(true).mixer_max_batch_size(),
        4_999_000
    );
    assert_eq!(
        SplitSizes::new(false, Some(1001), None).mixer_max_batch_size(),
        1
    );
    assert_eq!(
        SplitSizes::new(false, Some(500), None).mixer_max_batch_size(),
        -500
    );
    assert_eq!(
        SplitSizes::new(false, Some(-5), None).mixer_max_batch_size(),
        -1005
    );
    assert_eq!(
        SplitSizes::new(false, Some(i32::MIN), None).mixer_max_batch_size(),
        i32::MAX - 999
    );
}

#[test]
fn small_chunks_cover_the_queries_with_the_overlap() {
    let sizes = SplitSizes::new(false, Some(1200), Some(50));
    let chunks = split_query_batch(&[1000, 1500, 700], sizes, QueryStrand::Both).expect("split");
    assert!(chunks.len() > 1);
    for chunk in &chunks {
        assert!(!chunk.queries.is_empty());
    }
}

// Oracle (E2g audits b and round 2): a negative chunk size above a negative overlap
// splits a batch only when the wrapped difference fits in it at least twice; -1 over
// INT_MIN does not split 12000 residues (NCBI searches as without the variables), -5
// over -10 does.
#[test]
fn a_negative_chunk_size_splits_only_a_long_enough_batch() {
    let wide = SplitSizes::new(false, Some(-1), Some(i32::MIN));
    assert!(wide.negative_chunk_size());
    assert_eq!(calculate_num_chunks(wide, 12_000).0, 1);
    let narrow = SplitSizes::new(false, Some(-5), Some(-10));
    assert!(calculate_num_chunks(narrow, 12_000).0 > 1);
    assert!(!SplitSizes::new(false, Some(3000), Some(-10)).negative_chunk_size());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:149,166-170
// The chunk ranges in 64-bit `size_t` without debug-build underflow or overflow:
// - positive pairs with evened chunk size below overlap: (140, 100) on 90 bases (evened 96),
//   (120, 100) on 40 bases, (101, 100) on 2 bases, and (100, 95) on 31 bases (evened 85) wrap
//   `chunk_start` past the batch after chunk 0, leaving chunks 1.. empty;
// - positive pair at the exact overlap boundary (150, 100): 99 bases does not split, 100 bases
//   gives evened chunk size 100 and stops after chunk 0 leaving chunk 1 empty, 101 bases gives
//   two non-empty overlapping chunks (0..100 and 0..101, NCBI re-split -> Exception 1);
// - negative pairs: (-100, -200) on 310 bases wraps `chunk_size` past 310 (chunks 1 and 2 empty),
//   (-133, -200) on 200 bases and (-132, -200) on 201 bases keep `chunk_size` 1 without wrapping
//   yet advance `chunk_start` past or to the batch end (chunk 1 empty), while (-100, -200) on 248
//   bases and (-1, -200) on 930 bases leave 200-base gaps with every chunk holding a query.
#[test]
fn the_chunk_ranges_wrap_as_ncbis_size_t() {
    let ranges = |lengths: &[usize], sizes| -> Vec<Vec<(usize, usize)>> {
        split_query_batch(lengths, sizes, QueryStrand::Both)
            .unwrap()
            .iter()
            .map(|chunk| {
                chunk
                    .queries
                    .iter()
                    .map(|part| (part.from, part.to))
                    .collect()
            })
            .collect()
    };
    let negative = SplitSizes::new(false, Some(-100), Some(-200));
    assert_eq!(
        ranges(&[310], negative),
        vec![vec![(0, 310)], vec![], vec![]]
    );
    assert_eq!(
        ranges(&[248], negative),
        vec![vec![(0, 25)], vec![(225, 248)]]
    );
    assert_eq!(
        ranges(&[200], SplitSizes::new(false, Some(-133), Some(-200))),
        vec![vec![(0, 1)], vec![]]
    );
    assert_eq!(
        ranges(&[201], SplitSizes::new(false, Some(-132), Some(-200))),
        vec![vec![(0, 1)], vec![]]
    );
    assert_eq!(
        ranges(&[930], SplitSizes::new(false, Some(-1), Some(-200))),
        vec![
            vec![(0, 83)],
            vec![(283, 366)],
            vec![(566, 649)],
            vec![(849, 930)],
        ]
    );
    let close = SplitSizes::new(false, Some(100), Some(95));
    assert_eq!(
        ranges(&[31], close),
        vec![vec![(0, 31)], vec![], vec![], vec![], vec![], vec![]]
    );
    assert_eq!(
        ranges(&[90], SplitSizes::new(false, Some(140), None)),
        vec![vec![(0, 90)], vec![]]
    );
    assert_eq!(
        ranges(&[40], SplitSizes::new(false, Some(120), None)),
        vec![vec![(0, 40)], vec![]]
    );
    assert_eq!(
        ranges(&[2], SplitSizes::new(false, Some(101), None)),
        vec![vec![(0, 2)], vec![]]
    );
    let c150 = SplitSizes::new(false, Some(150), None);
    assert!(split_query_batch(&[99], c150, QueryStrand::Both).is_none());
    assert_eq!(ranges(&[100], c150), vec![vec![(0, 100)], vec![]]);
    assert_eq!(ranges(&[101], c150), vec![vec![(0, 100)], vec![(0, 101)]]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:392-410,604-626
// A strand that is not searched has no chunk context and the offset INT4_MAX; the other
// strand's offsets are those of the split with both strands.
#[test]
fn a_single_strand_split_keeps_the_offsets_of_the_strand_it_searches() {
    let lengths = [1_999_976, 1_000_000];
    let sizes = SplitSizes::default_for(false);
    let both = split_query_batch(&lengths, sizes, QueryStrand::Both).expect("split");
    for (strand, keep) in [(QueryStrand::Plus, 0usize), (QueryStrand::Minus, 1)] {
        let one = split_query_batch(&lengths, sizes, strand).expect("split");
        assert_eq!(one.len(), both.len());
        for (one, both) in one.iter().zip(&both) {
            assert_eq!(one.queries, both.queries);
            assert_eq!(one.contexts.len(), both.contexts.len());
            for (position, (&context, &offset)) in
                one.contexts.iter().zip(&one.context_offsets).enumerate()
            {
                if position % 2 == keep {
                    assert_eq!(context, both.contexts[position]);
                    assert_eq!(offset, both.context_offsets[position]);
                } else {
                    assert_eq!(context, None);
                    assert_eq!(offset, i32::MAX);
                }
            }
        }
    }
    let plus = split_query_batch(&lengths, sizes, QueryStrand::Plus).expect("split");
    assert_eq!(plus[1].contexts, vec![Some(0), None, Some(2), None]);
    let minus = split_query_batch(&lengths, sizes, QueryStrand::Minus).expect("split");
    assert_eq!(minus[1].contexts, vec![None, Some(1), None, Some(3)]);
    assert_eq!(
        minus[1].context_offsets,
        vec![i32::MAX, 0, i32::MAX, 1_000_000 - 42]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:233-290
// `SplitBatchSteps` with `QueryChunk` and `without_resplit_stop()` (BLASTN): empty chunks refuse
// through `empty_chunk_refusal`, re-splitting chunks (Exception 1: 150/100 on L=101) proceed,
// `ReadFails::Setup` stops before validity, `setup_checks` run once on the first valid chunk (or
// become `UnmetChecks` when no chunk is valid), `ReadFails::Search` stops after `setup_checks`,
// and D13 negative-pair gaps refuse only with subjects in hand.
#[test]
fn blastn_split_batch_steps_follow_ncbis_order_and_skip_resplit_stop() {
    let stop = |error: anyhow::Error| error;
    assert_eq!(
        empty_chunk_refusal("BLASTN"),
        "a query chunk without a query is not supported by LOSAT's BLASTN"
    );
    // Empty chunk (140/100 on L=90) refuses in `new` before any read or check.
    let empty_sizes = SplitSizes::new(false, Some(140), None);
    let empty_chunks = split_query_batch(&[90], empty_sizes, QueryStrand::Both).expect("split");
    let err = SplitBatchSteps::new(
        "BLASTN",
        &empty_chunks,
        empty_sizes,
        Some(ReadFails::Setup),
        Err(SetupCheckFailure::Engine(anyhow::anyhow!("short"))),
        &stop,
    )
    .err()
    .expect("empty chunk stops");
    assert_eq!(err.to_string(), empty_chunk_refusal("BLASTN"));

    // Approved Exception 1 (150/100 on L=101): `chunk_would_be_split` is true, so default stops
    // while `without_resplit_stop` searches both chunks.
    let ex1_sizes = SplitSizes::new(false, Some(150), None);
    let ex1_chunks = split_query_batch(&[101], ex1_sizes, QueryStrand::Both).expect("split");
    assert!(chunk_would_be_split(&ex1_chunks, ex1_sizes));
    let mut strict = SplitBatchSteps::new("BLASTN", &ex1_chunks, ex1_sizes, None, Ok(()), &stop)
        .expect("non-empty");
    assert!(strict.before_read(&ex1_chunks[0], true).is_err());
    let mut blastn = SplitBatchSteps::new("BLASTN", &ex1_chunks, ex1_sizes, None, Ok(()), &stop)
        .expect("non-empty")
        .without_resplit_stop();
    for chunk in &ex1_chunks {
        assert_eq!(blastn.before_read(chunk, true).unwrap(), ChunkStep::Read);
        blastn.before_search(true).expect("positive pair");
    }

    // Negative pair with gaps (-100/-200 on L=248): OID 0 read stops before validity; invalid
    // chunks leave `unmet_checks` (`TracebackFailed` for engine error, preserved refusal for
    // count-0/greedy); valid chunk meets setup checks, then `ReadFails::Search`, then D13 gaps.
    let gap_sizes = SplitSizes::new(false, Some(-100), Some(-200));
    let gap_chunks = split_query_batch(&[248], gap_sizes, QueryStrand::Both).expect("split");
    let mut setup_read = SplitBatchSteps::new(
        "BLASTN",
        &gap_chunks,
        gap_sizes,
        Some(ReadFails::Setup),
        Err(SetupCheckFailure::Engine(anyhow::anyhow!("short"))),
        &stop,
    )
    .unwrap()
    .without_resplit_stop();
    assert_eq!(
        setup_read.before_read(&gap_chunks[0], false).unwrap(),
        ChunkStep::End
    );

    let mut gapvalid_engine = SplitBatchSteps::new(
        "BLASTN",
        &gap_chunks,
        gap_sizes,
        Some(ReadFails::Search),
        Err(SetupCheckFailure::Engine(anyhow::anyhow!("short"))),
        &stop,
    )
    .unwrap()
    .without_resplit_stop();
    for chunk in &gap_chunks {
        assert_eq!(
            gapvalid_engine.before_read(chunk, false).unwrap(),
            ChunkStep::Skip
        );
    }
    assert!(matches!(
        gapvalid_engine.unmet_checks(),
        Some(UnmetChecks::TracebackFailed)
    ));

    let mut gapvalid_refusal = SplitBatchSteps::new(
        "BLASTN",
        &gap_chunks,
        gap_sizes,
        None,
        Err(SetupCheckFailure::Refusal(anyhow::anyhow!("count 0"))),
        &stop,
    )
    .unwrap()
    .without_resplit_stop();
    for chunk in &gap_chunks {
        assert_eq!(
            gapvalid_refusal.before_read(chunk, false).unwrap(),
            ChunkStep::Skip
        );
    }
    assert!(matches!(
        gapvalid_refusal.unmet_checks(),
        Some(UnmetChecks::Refusal(error)) if error.to_string() == "count 0"
    ));

    let mut gap_valid = SplitBatchSteps::new("BLASTN", &gap_chunks, gap_sizes, None, Ok(()), &stop)
        .unwrap()
        .without_resplit_stop();
    assert_eq!(
        gap_valid.before_read(&gap_chunks[0], true).unwrap(),
        ChunkStep::Read
    );
    gap_valid.before_search(false).expect("no subjects in hand");
    let gap_err = gap_valid.before_search(true).unwrap_err();
    assert_eq!(gap_err.to_string(), negative_split_refusal("BLASTN"));
}
