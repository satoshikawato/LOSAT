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
