use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:123-138
// A 60000-residue tblastn query is 3 chunks of 20067 residues that overlap by 100.
#[test]
fn a_long_tblastn_query_is_three_overlapping_chunks() {
    let sizes = protein_split_sizes(20000, None, None);
    let chunks = split_protein_batch(&[60_000], sizes).unwrap();
    let ranges: Vec<_> = chunks
        .iter()
        .map(|chunk| (chunk.queries[0].from, chunk.queries[0].to))
        .collect();
    assert_eq!(ranges, vec![(0, 20067), (19967, 40034), (39934, 60000)]);
    let offsets: Vec<_> = chunks
        .iter()
        .map(|chunk| chunk.context_offsets[0])
        .collect();
    assert_eq!(offsets, vec![0, 19967, 39934]);
    assert!(split_protein_batch(&[39_799], sizes).is_none());
    assert!(!protein_chunk_would_be_split(&chunks, sizes));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:149,166-170
// The chunk ranges in `size_t`, where a debug build's `usize` subtraction panicked:
// CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200 make 310 residues 3 chunks of 2^64 / 3 or so, the
// first of which covers the batch (and would be split again, but NCBI's `CQuerySplitter::Split`
// stops first at the chunks without a query) while the others keep NCBI's empty bounds, and 248
// residues 2 chunks of 25 with a gap of 200 between them; CHUNK_SIZE=100
// OVERLAP_CHUNK_SIZE=95 make 31 residues 6 chunks of 85, below the overlap, the first covering
// the batch.
#[test]
fn the_chunk_ranges_wrap_as_ncbis_size_t() {
    let ranges = |lengths: &[usize], sizes| -> Vec<Vec<(usize, usize)>> {
        split_protein_batch(lengths, sizes)
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
    let negative = protein_split_sizes(20000, Some(-100), Some(-200));
    assert_eq!(
        ranges(&[310], negative),
        vec![vec![(0, 310)], vec![], vec![]]
    );
    let chunks = split_protein_batch(&[310], negative).unwrap();
    assert!(protein_chunk_would_be_split(&chunks, negative));
    assert_eq!(
        ranges(&[248], negative),
        vec![vec![(0, 25)], vec![(225, 248)]]
    );
    let chunks = split_protein_batch(&[248], negative).unwrap();
    assert!(!protein_chunk_would_be_split(&chunks, negative));
    let close = protein_split_sizes(20000, Some(100), Some(95));
    assert_eq!(
        ranges(&[31], close),
        vec![vec![(0, 31)], vec![], vec![], vec![], vec![], vec![]]
    );
}

// Two queries share the first chunk (10051 residues, `calculate_num_chunks`); the
// second query's later part starts where its first part's overlap begins.
#[test]
fn a_query_over_two_chunks_has_the_offset_of_its_first_part() {
    let sizes = protein_split_sizes(10000, None, None);
    let chunks = split_protein_batch(&[5_000, 15_000], sizes).unwrap();
    assert_eq!(chunks.len(), 2);
    assert_eq!(
        chunks[0].queries,
        vec![
            ChunkQuery {
                query: 0,
                from: 0,
                to: 5_000
            },
            ChunkQuery {
                query: 1,
                from: 0,
                to: 5_051
            },
        ]
    );
    assert_eq!(chunks[0].context_offsets, vec![0, 0]);
    assert_eq!(
        chunks[1].queries,
        vec![ChunkQuery {
            query: 1,
            from: 4_951,
            to: 15_000
        }]
    );
    assert_eq!(chunks[1].context_offsets, vec![4_951]);
}

// The lower-case mask [3, 6) restricted to the part [4, 8) is [4, 7) of the query (one
// residue longer); the mask [10, 12) is outside the part.
#[test]
fn a_part_keeps_ncbis_restricted_lower_case_masks() {
    let query = b"ACDefgHIKLmn";
    let part = ChunkQuery {
        query: 0,
        from: 4,
        to: 8,
    };
    assert_eq!(
        protein_chunk_part_in_record(query, 0, &part),
        b"fghI".to_vec()
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:278-281 (cited at
// `protein_chunk_part_in_record`): with `-query_loc 5-` of `MKVa...`, the run [3, 3] of the
// record grows to [3, 4] and masks the first searched letter in the first chunk; a run ending
// two letters before the interval does not reach it.
#[test]
fn a_run_before_the_interval_reaches_its_first_letter() {
    let part = ChunkQuery {
        query: 0,
        from: 0,
        to: 6,
    };
    assert_eq!(
        protein_chunk_part_in_record(b"MKVaWCDEFG", 4, &part),
        b"wCDEFG".to_vec()
    );
    assert_eq!(
        protein_chunk_part_in_record(b"MKaVWCDEFG", 4, &part),
        b"WCDEFG".to_vec()
    );
    assert_eq!(
        protein_chunk_part_in_record(b"MKVAWCDEFG", 4, &part),
        b"WCDEFG".to_vec()
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:233-290 (cited at
// `SplitBatchSteps`)
// The steps of a split batch before each chunk's search reads the subjects, in NCBI's order, as
// BLASTP's and TBLASTN's searches, the set-up of a refused batch and a call's setup failures meet
// them (review 4 of line RG item X-short-bat, M1; item X-bpgap): a chunk without a query before any
// chunk; per chunk a split again, then the setup's read of OID 0 (before the validity), then no
// valid context (the chunk is passed over), then the batch's setup checks at the first valid chunk
// (once), then the search's read; D13's gaps where a chunk's search would read subjects in hand;
// the checks that no chunk met after the last chunk.
#[test]
fn a_split_batch_meets_its_chunk_steps_in_ncbis_order() {
    #[derive(Debug)]
    struct Stop;
    impl std::fmt::Display for Stop {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str("stop")
        }
    }
    impl std::error::Error for Stop {}
    let stop = |refusal: anyhow::Error| refusal.context(Stop);
    let sizes = |chunk: i64, overlap: i64| SplitSizes {
        chunk_size: chunk as u64,
        overlap: overlap as u64,
    };
    let fails = || -> Result<(), SetupCheckFailure> {
        Err(SetupCheckFailure::Engine(anyhow::anyhow!("setup")))
    };
    let is_stop = |error: anyhow::Error| error.is::<Stop>();
    let is_setup_error = |error: anyhow::Error| !error.is::<Stop>() && error.to_string() == "setup";
    // `CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200`: a first chunk that covers a query of 310 residues,
    // then chunks without a query, before any read or check.
    let empty = split_protein_batch(&[310], sizes(-100, -200)).expect("split");
    let refused = SplitBatchSteps::new(
        "TBLASTN",
        &empty,
        sizes(-100, -200),
        Some(ReadFails::Setup),
        fails(),
        &stop,
    );
    assert!(is_stop(refused.err().expect("a chunk without a query")));
    // q1, q2, q1 of 300 residues with 300 and 10: 3 chunks of 307, none split again.
    let chunks = split_protein_batch(&[300, 300, 300], sizes(300, 10)).expect("split");
    assert_eq!(chunks.len(), 3);
    let steps = |read_fails, checks| {
        SplitBatchSteps::new("BLASTP", &chunks, sizes(300, 10), read_fails, checks, &stop)
            .expect("no empty chunk")
    };
    let mut failing = steps(None, fails());
    assert_eq!(
        failing.before_read(&chunks[0], false).unwrap(),
        ChunkStep::Skip
    );
    assert!(is_setup_error(
        failing.before_read(&chunks[1], true).unwrap_err()
    ));
    assert!(failing.unmet_checks().is_none());
    let mut passing = steps(None, Ok(()));
    for chunk in &chunks {
        assert_eq!(passing.before_read(chunk, true).unwrap(), ChunkStep::Read);
        passing
            .before_search(true)
            .expect("a positive pair has no gap");
    }
    // The setup's read of OID 0 ends the stage before the chunk's validity and the checks.
    let mut setup_read = steps(Some(ReadFails::Setup), fails());
    assert_eq!(
        setup_read.before_read(&chunks[0], false).unwrap(),
        ChunkStep::End
    );
    // The search's read comes after the checks of the first valid chunk.
    let mut search_read = steps(Some(ReadFails::Search), fails());
    assert_eq!(
        search_read.before_read(&chunks[0], false).unwrap(),
        ChunkStep::Skip
    );
    assert!(is_setup_error(
        search_read.before_read(&chunks[1], true).unwrap_err()
    ));
    let mut search_read = steps(Some(ReadFails::Search), Ok(()));
    assert_eq!(
        search_read.before_read(&chunks[0], true).unwrap(),
        ChunkStep::End
    );
    // No chunk with a valid context: nothing is read, and the traceback's own setup meets the
    // checks that no chunk met.
    let mut unread = steps(Some(ReadFails::Search), fails());
    assert!(!split_batch_reaches_read(&mut unread, &chunks, |_| Ok(false)).unwrap());
    assert!(matches!(
        unread.unmet_checks(),
        Some(UnmetChecks::TracebackFailed)
    ));
    // LOSAT's refusal there (a count of 0) stays the refusal, at the first valid chunk and after
    // the last chunk.
    let refusal = || -> Result<(), SetupCheckFailure> {
        Err(SetupCheckFailure::Refusal(anyhow::anyhow!("count 0")))
    };
    let mut refused = steps(None, refusal());
    let error = refused.before_read(&chunks[0], true).unwrap_err();
    assert!(!error.is::<Stop>() && error.to_string() == "count 0");
    let mut refused = steps(None, refusal());
    assert!(!split_batch_reaches_read(&mut refused, &chunks, |_| Ok(false)).unwrap());
    assert!(matches!(
        refused.unmet_checks(),
        Some(UnmetChecks::Refusal(error)) if error.to_string() == "count 0"
    ));
    let mut unread = steps(Some(ReadFails::Search), Ok(()));
    assert!(!split_batch_reaches_read(&mut unread, &chunks, |_| Ok(false)).unwrap());
    assert!(unread.unmet_checks().is_none());
    let mut read = steps(Some(ReadFails::Search), Ok(()));
    assert!(split_batch_reaches_read(&mut read, &chunks, |chunk| Ok(chunk == &chunks[2])).unwrap());
    // 453 and 227: chunk 0 of the same batch would be split again, which comes before the setup's
    // read and the checks (the set-up records q1 and q2 alone make 2 chunks that are not).
    let again = split_protein_batch(&[300, 300, 300], sizes(453, 227)).expect("split");
    let mut again_steps = SplitBatchSteps::new(
        "TBLASTN",
        &again,
        sizes(453, 227),
        Some(ReadFails::Setup),
        fails(),
        &stop,
    )
    .expect("no empty chunk");
    assert!(is_stop(
        again_steps.before_read(&again[0], true).unwrap_err()
    ));
    let pair = split_protein_batch(&[300, 300], sizes(453, 227)).expect("split");
    assert!(!protein_chunk_would_be_split(&pair, sizes(453, 227)));
    // `CHUNK_SIZE=-1 OVERLAP_CHUNK_SIZE=-200` on 398 residues: chunks [0, 100) and [300, 398),
    // each with a query, with the gap between them. D13 refuses where a chunk's search would read
    // subjects in hand, and without one the batch is searched on.
    let gaps = split_protein_batch(&[398], sizes(-1, -200)).expect("split");
    let ranges: Vec<_> = gaps
        .iter()
        .map(|chunk| (chunk.queries[0].from, chunk.queries[0].to))
        .collect();
    assert_eq!(ranges, vec![(0, 100), (300, 398)]);
    let mut gap_steps = SplitBatchSteps::new("BLASTP", &gaps, sizes(-1, -200), None, Ok(()), &stop)
        .expect("each chunk holds a query");
    assert!(!gap_steps.gap_step_pending(false));
    assert!(gap_steps.gap_step_pending(true));
    assert_eq!(
        gap_steps.before_read(&gaps[0], true).unwrap(),
        ChunkStep::Read
    );
    gap_steps
        .before_search(false)
        .expect("no subject in hand: nothing is lost");
    assert!(is_stop(gap_steps.before_search(true).unwrap_err()));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:60-61 (cited at
// `splits_query_batches`)
// The one rule of whether NCBI splits a search's query batches (review 1 of line RG item X-bpgap,
// S6): a gapped search of a FASTA query against sequences only, without a PHI-BLAST pattern
// (psiblast `-phi_pattern`: NCBI 2.17.0 searches a 24,334-residue query whole, line PS item PS8a
// `PS/tmp/fix8a/lq/ncbi.q_24k.6`).
#[test]
fn ncbi_splits_only_a_gapped_search_of_a_fasta_query_against_sequences() {
    assert!(splits_query_batches(true, false, false, false));
    assert!(!splits_query_batches(false, false, false, false));
    assert!(!splits_query_batches(true, true, false, false));
    assert!(!splits_query_batches(true, false, true, false));
    assert!(!splits_query_batches(true, false, false, true));
}
