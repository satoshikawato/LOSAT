use super::*;
use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::config::ScoringMatrix;
use crate::stats::lookup_protein_params_ungapped;
use crate::utils::matrix::ncbistdaa;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:87-129
fn make_poly_a_query_frame() -> QueryFrame {
    QueryFrame {
        frame: 1,
        aa_seq: vec![
            0,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            0,
        ],
        aa_seq_nomask: None,
        aa_len: 6,
        orig_len: 6,
        seg_masks: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-456
fn make_poly_a_lookup() -> BlastAaLookupTable {
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let queries = vec![vec![make_poly_a_query_frame()]];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2746-2748
    let (lookup, _) = build_ncbi_lookup(&queries, 0, &ungapped, true);
    lookup
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aascan_unit_test.cpp:288-326 ScanOffsetTest
#[test]
fn test_s_blast_aa_scan_subject_small_buffer_resumes_without_splitting_subject_offset() {
    let lookup = make_poly_a_lookup();
    let subject = vec![ncbistdaa::A; 6];
    let seq_ranges = [(0, subject.len() as i32)];
    let mut offset_pairs = vec![OffsetPair::default(); 5];
    let array_size = i32::try_from(offset_pairs.len())
        .expect("NCBI BLAST offset pair array size must fit in Int4");
    let mut scan_range = [0, 0, subject.len() as i32 - lookup.word_length];

    let mut observed_subject_offsets = Vec::new();
    while scan_range[1] <= scan_range[2] {
        let hits = s_blast_aa_scan_subject(
            &lookup,
            &subject,
            &seq_ranges,
            &mut offset_pairs,
            array_size,
            &mut scan_range,
        );
        assert_eq!(hits, 4);

        let batch = &offset_pairs[..hits as usize];
        let expected_subject_offset = observed_subject_offsets.len() as u32;
        assert!(batch
            .iter()
            .all(|pair| pair.s_off == expected_subject_offset));
        assert_eq!(
            batch.iter().map(|pair| pair.q_off).collect::<Vec<_>>(),
            vec![0, 1, 2, 3]
        );
        observed_subject_offsets.push(expected_subject_offset);
    }

    assert_eq!(observed_subject_offsets, vec![0, 1, 2, 3]);
}
