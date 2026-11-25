//! Focused Stage C boundary checks, source-derived expectations.
use super::{
    args::ResolvedOptions,
    input::{read_fasta, FastaRecord},
    parameters::{effective_lengths, score_block},
    query_setup::{prepare_queries, PreparedQueryBatch},
    split::{
        calculate_num_chunks, chunk_stops, out_of_frame, split_queries, validate_chunk_size,
        QuerySplitRefusal, DEFAULT_SIZES,
    },
};
use crate::algorithm::blastn::query_split::SplitSizes;
use crate::cli::Commands;
use std::path::Path;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:49-62
pub(super) fn options() -> ResolvedOptions {
    let cli: crate::cli::Cli = crate::cli::try_parse_from([
        "losat",
        "blastx",
        "-query",
        "unused.fna",
        "-subject",
        "unused.faa",
        "-seg",
        "no",
        "-comp_based_stats",
        "0",
    ])
    .unwrap();
    let Commands::Blastx(args) = cli.command else {
        unreachable!()
    };
    args.resolve().unwrap()
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-264
pub(super) fn fixture(name: &str) -> (Vec<FastaRecord>, PreparedQueryBatch, ResolvedOptions) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../docs/evidence/losatx_stage_a/inputs/generated")
        .join(name);
    let records = read_fasta(&path, false, false).unwrap();
    let o = options();
    let batch = prepare_queries(&records, &o).unwrap();
    (records, batch, o)
}
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/split_query_unit_test.cpp:1901-1925 CalculateNumberChunks
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/split_query_unit_test.cpp:1927-1931 InvalidChunkSizeBlastx
#[test]
fn active_n04_mutable_chunk_size_and_invalid_divisibility() {
    let mut size: u64 = 10002;
    assert_eq!(calculate_num_chunks(&mut size, 297, 10240000, 1), 1055);
    let length = size as usize / 2;
    assert_eq!(calculate_num_chunks(&mut size, 297, length, 1), 1);
    let length = 3 * size as usize - 2 * 297;
    assert_eq!(calculate_num_chunks(&mut size, 297, length, 1), 3);
    let length = 1 + 2 * size as usize + 297;
    assert_eq!(calculate_num_chunks(&mut size, 297, length, 1), 2);
    assert!(validate_chunk_size(40000).is_err());
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:748-791
#[test]
fn split_modulo_frame_maps_and_shared_positive_corrections() {
    for (len, map, negative) in [
        (19409, [0, 1, 2, 3, 4, 5], [0, 0, 0]),
        (19410, [0, 1, 2, 3, 4, 5], [3136, 3136, 3136]),
        (19411, [0, 1, 2, 4, 5, 3], [3136, 3136, 3137]),
    ] {
        let (records, full, o) = fixture(&format!("S12_split{len}.fna"));
        let chunks = split_queries(&records, &full, &o, DEFAULT_SIZES).unwrap();
        assert_eq!(chunks.len(), if len < 19410 { 1 } else { 2 });
        assert_eq!(chunks[0].absolute_contexts, map);
        assert_eq!(&chunks[0].corrections[3..], negative);
        if chunks.len() > 1 {
            assert_eq!(chunks[1].range, 9705..len);
            assert_eq!(chunks[1].corrections, [3235, 3235, 3235, 0, 0, 0]);
        }
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:76-100
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:56-61
#[test]
fn multiple_queries_and_ungapped_do_not_split() {
    let (mut records, full, mut o) = fixture("S12_split19411.fna");
    o.gapped = false;
    assert_eq!(
        split_queries(&records, &full, &o, DEFAULT_SIZES)
            .unwrap()
            .len(),
        1
    );
    o.gapped = true;
    records.push(records[0].clone());
    let full = prepare_queries(&records, &o).unwrap();
    assert_eq!(
        split_queries(&records, &full, &o, DEFAULT_SIZES)
            .unwrap()
            .len(),
        1
    );
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:145-171
// The sizes of CHUNK_SIZE and OVERLAP_CHUNK_SIZE (as `split_sizes` reads them) set the chunks;
// the sizes NCBI cannot split with are LOSAT's refusals.
#[test]
fn chunk_ranges_follow_the_sizes_and_unsplittable_sizes_are_refused() {
    let (records, full, o) = fixture("S12_split19411.fna");
    let split = |chunk_size, overlap| {
        split_queries(
            &records,
            &full,
            &o,
            SplitSizes {
                chunk_size,
                overlap,
            },
        )
    };
    let chunks = split(3000, 0).unwrap();
    let ranges: Vec<_> = chunks.iter().map(|chunk| chunk.range.clone()).collect();
    assert_eq!(
        ranges,
        [
            0..3000,
            3000..6000,
            6000..9000,
            9000..12000,
            12000..15000,
            15000..19411
        ]
    );
    assert!(chunks.iter().all(|chunk| chunk.overlap == 0));
    let chunks = split(3000, 297).unwrap();
    assert_eq!(chunks.len(), 7);
    assert_eq!(chunks[1].range, 2703..5703);
    assert_eq!(chunks[6].range, 16218..19411);
    assert!(chunks.iter().all(|chunk| chunk.overlap == 99));
    // Not above the overlap (CHUNK_SIZE 0, 3, -1 with the default overlap; a negative
    // overlap), or one chunk of the whole query.
    for (chunk_size, overlap) in [
        (0, 297),
        (3, 297),
        (u64::MAX, 297),
        (10002, u64::MAX - 2),
        (30000, 297),
    ] {
        assert_eq!(split(chunk_size, overlap).unwrap().len(), 1);
    }
    // 2 x the overlap at least the chunk size; -1 above -4 (chunks without a query): the stops
    // before the chunks are set up, also through `chunk_stops`.
    for (chunk_size, overlap) in [(600, 300), (10002, 5001), (u64::MAX, u64::MAX - 3)] {
        let sizes = SplitSizes {
            chunk_size,
            overlap,
        };
        let error = split(chunk_size, overlap).unwrap_err();
        assert!(
            error.is::<QuerySplitRefusal>(),
            "{chunk_size} {overlap}: {error}"
        );
        let stop = chunk_stops(&records, &full, &o, sizes).unwrap_err();
        assert_eq!(stop.to_string(), error.to_string());
    }
    assert!(chunk_stops(&records, &full, &o, DEFAULT_SIZES).is_ok());
    // An overlap not divisible by 3 is split as NCBI splits it, and `out_of_frame` refuses the
    // chunks after the first; a multiple of 3 or one chunk is not refused.
    let sizes = |overlap| SplitSizes {
        chunk_size: 3000,
        overlap,
    };
    let chunks = split(3000, 100).unwrap();
    assert!(chunks.len() > 1);
    let error = out_of_frame(&chunks, sizes(100), &o).expect("refused");
    assert!(error.is::<QuerySplitRefusal>(), "{error}");
    assert!(out_of_frame(&split(3000, 99).unwrap(), sizes(99), &o).is_none());
    assert!(out_of_frame(&chunks[..1], sizes(100), &o).is_none());
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:776-785
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:846-847
#[test]
fn invalid_statistics_keep_fixed_space_with_zero_adjustment() {
    let (_, full, o) = fixture("S06_seed.fna");
    let mut p = score_block(&full);
    p[0].valid = false;
    let spaces = vec![37; full.contexts.len()];
    effective_lengths(&full, &mut p, &o, 83, 1, Some(&spaces)).unwrap();
    assert_eq!(p[0].search_space, 37);
    assert_eq!(p[0].length_adjustment, 0);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2781-2787
#[test]
fn statistical_validity_is_distinct_from_preparation_validity() {
    let (_, mut full, _) = fixture("S06_seed.fna");
    assert!(full.contexts.iter().all(|c| c.is_valid));
    full.sequence_start.fill(0);
    let p = score_block(&full);
    assert!(p.iter().all(|p| !p.valid));
    assert!(full.contexts.iter().all(|c| c.is_valid));
}
/// No NCBI counterpart (`LOSAT_X_BXKAONCE`): for a batch that is not split, the blocks `score_block`
/// gave before `effective_lengths` (what `search_core` keeps) equal `score_block` of the chunk
/// `split_queries` returns. Several queries with contexts that are invalid (a query of 4 bases,
/// one of N only: its contexts are prepared valid and statistically invalid) among valid ones,
/// gapped or not; and one query of 19,410 letters or more without gapped (`split_queries` makes
/// one chunk, `!options.gapped`).
#[test]
fn kept_blocks_equal_score_block_of_the_unsplit_chunk() {
    let mut state = 12345u64;
    let mut random = |letters: usize| -> Vec<u8> {
        (0..letters)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 60) as usize & 3]
            })
            .collect()
    };
    let record = |sequence: Vec<u8>| FastaRecord {
        internal_id: String::new(),
        title: Vec::new(),
        sequence,
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    };
    let sequence = random(600);
    let several = vec![
        record(sequence.clone()),
        record(b"ACGT".to_vec()),
        record(vec![b'N'; 300]),
        record(sequence[100..400].to_vec()),
    ];
    let long = vec![record(random(19_500))];
    for (records, gapped) in [(&several, true), (&several, false), (&long, false)] {
        let mut o = options();
        o.gapped = gapped;
        let mut full = prepare_queries(records, &o).unwrap();
        let mut params = score_block(&full);
        if records.len() > 1 {
            assert!(params.iter().any(|p| p.valid) && params.iter().any(|p| !p.valid));
            // The assignment below flips a context's `is_valid` from true to false.
            assert!(full
                .contexts
                .iter()
                .zip(&params)
                .any(|(c, p)| c.is_valid && !p.valid));
        }
        let kept = params.clone();
        effective_lengths(&full, &mut params, &o, 1_000_000, 100, None).unwrap();
        for (context, p) in full.contexts.iter_mut().zip(&params) {
            context.is_valid = p.valid;
        }
        let chunks = split_queries(records, &full, &o, DEFAULT_SIZES).unwrap();
        assert_eq!(chunks.len(), 1);
        let again = score_block(&chunks[0].prepared);
        assert_eq!(kept.len(), again.len());
        for (a, b) in kept.iter().zip(&again) {
            assert_eq!(a.valid, b.valid);
            for (x, y) in [
                (a.ungapped.lambda, b.ungapped.lambda),
                (a.ungapped.k, b.ungapped.k),
                (a.ungapped.h, b.ungapped.h),
                (a.ungapped.alpha, b.ungapped.alpha),
                (a.ungapped.beta, b.ungapped.beta),
            ] {
                assert_eq!(x.to_bits(), y.to_bits());
            }
            assert_eq!(a.search_space, b.search_space);
            assert_eq!(a.length_adjustment, b.length_adjustment);
            assert_eq!(a.word_cutoff, b.word_cutoff);
            assert_eq!(a.word_xdrop_init, b.word_xdrop_init);
            assert_eq!(a.word_xdrop, b.word_xdrop);
            assert_eq!(a.hit_cutoff, b.hit_cutoff);
            assert_eq!(a.hit_cutoff_max, b.hit_cutoff_max);
        }
    }
}
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:282-292 BackboneSequenceTest
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:294-304 SmallboneSequenceTest
#[test]
fn active_ncbi_lookup_overflow_u16_boundary_preserves_all_offsets() {
    use crate::algorithm::tblastx::lookup::{build_ncbi_lookup_from_prepared, QueryContext};
    for len in [65533usize, 65534] {
        let c = QueryContext {
            q_idx: 0,
            frame: 1,
            aa_seq: vec![0; len + 2],
            aa_seq_nomask: None,
            aa_len: len,
            orig_len: len * 3,
            frame_base: 0,
            is_valid: true,
            karlin_params: Default::default(),
        };
        let lookup = build_ncbi_lookup_from_prepared(
            vec![0; len + 1],
            vec![(0, len as i32 - 1)],
            vec![c],
            0,
            false,
            1,
            None,
        );
        let crate::algorithm::tblastx::lookup::XBxTable::Wide(lookup) = lookup else {
            panic!("a build without a cutoff takes Int4 hits");
        };
        let hits = lookup.get_hits(0);
        assert_eq!(hits.len(), len - 2);
        assert_eq!(hits.last(), Some(&(len as i32 - 3)));
        assert!(hits.iter().copied().eq(0..len as i32 - 2));
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250
#[test]
fn unported_long_subject_chunking_is_rejected_before_search() {
    let (records, _, o) = fixture("S06_seed.fna");
    let subject = FastaRecord {
        internal_id: "Subject_1".into(),
        title: b"long".to_vec(),
        sequence: vec![b'W'; 5_000_001],
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    };
    let error = super::search::search_preliminary(&records, &[subject], &o)
        .expect_err("unported subject splitting");
    assert_eq!(
        error.to_string(),
        "subject sequence 1 has 5000001 residues; a subject longer than 5000000 residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTX"
    );
}

/// A batch without a valid Karlin-Altschul context is not searched, so it does not reach the
/// subjects: NCBI writes its report, and LOSAT does not refuse a long subject for it. The subject is
/// the one that `unported_long_subject_chunking_is_rejected_before_search` refuses when its batch is
/// searched, so `Ok` shows that this batch was not.
#[test]
fn a_batch_that_is_not_searched_does_not_refuse_a_long_subject() {
    let (_, _, o) = fixture("S06_seed.fna");
    let record = |id: &str, sequence: Vec<u8>| FastaRecord {
        internal_id: id.into(),
        title: Vec::new(),
        sequence,
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    };
    let query = record("Query_1", vec![b'N'; 300]);
    let subject = record("Subject_1", vec![b'W'; 5_000_001]);
    super::search::search_preliminary(&[query], &[subject], &o).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:233-239
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205-205
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-987
#[test]
fn short_subject_cutoff_uses_adapter_minimum_ten() {
    let o = options();
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../docs/evidence/losatx_stage_c/inputs/compressed_overflow.fna");
    let records = read_fasta(&path, false, false).unwrap();
    let full = prepare_queries(&records, &o).unwrap();
    let mut p = score_block(&full);
    super::parameters::subject_parameters(
        &full,
        &mut p,
        &o,
        super::parameters::SubjectNumbers::one(8),
        None,
    )
    .unwrap();
    // Pinned NCBI CALL0 context0: qlen6468, subject8, min length10 -> max cutoff18.
    assert_eq!(p[0].hit_cutoff_max, 18);
    assert_eq!(p[0].word_cutoff, 18);
}
