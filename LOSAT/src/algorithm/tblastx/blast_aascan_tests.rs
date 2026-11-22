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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58,
// c++/src/algo/blast/core/blast_aascan.c:75-123 and c++/src/algo/blast/api/seqsrc_seqdb.cpp:398-405
// The ranges of masks stored out of order, `[18, 21)` before `[6, 9)`: `[0, 18)`, `[21, 6)` (no
// word) and `[9, 30)`. The scan, resumed after each word of four hits, gives the words 0 to 15 and
// then 9 to 27: the subject offset falls, and `ScanOrder` takes every hit in that order.
#[test]
fn scan_order_takes_the_hits_of_ranges_out_of_order_as_the_scan_gives_them() {
    let lookup = make_poly_a_lookup();
    let subject = vec![ncbistdaa::A; 30];
    let seq_ranges = [(0, 18), (21, 6), (9, 30)];
    let mut offset_pairs = vec![OffsetPair::default(); 5];
    let array_size = i32::try_from(offset_pairs.len())
        .expect("NCBI BLAST offset pair array size must fit in Int4");
    let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - lookup.word_length];

    let mut subject_offsets = Vec::new();
    while scan_range[1] <= scan_range[2] {
        let hits = s_blast_aa_scan_subject(
            &lookup,
            &subject,
            &seq_ranges,
            &mut offset_pairs,
            array_size,
            &mut scan_range,
        );
        subject_offsets.extend(offset_pairs[..hits as usize].iter().map(|pair| pair.s_off));
    }

    let words: Vec<u32> = (0..=15).chain(9..=27).collect();
    let expected: Vec<u32> = words.iter().flat_map(|&s_off| [s_off; 4]).collect();
    assert_eq!(subject_offsets, expected);
    assert!(subject_offsets.windows(2).any(|pair| pair[1] < pair[0]));
    let mut order = ScanOrder::default();
    assert!(subject_offsets.iter().all(|&s_off| order.takes(
        &seq_ranges,
        lookup.word_length,
        s_off as i32
    )));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58 and
// c++/src/algo/blast/core/aa_ungapped.c:497-500
// A hit before the last one of its range, past the range's last word or in a range without a
// word is not in the scan's order; the first range, shorter than a word, has the word at its left.
#[test]
fn scan_order_refuses_hits_the_scan_cannot_give_next() {
    let ranges = [(0, 30)];
    let mut order = ScanOrder::default();
    assert!(order.takes(&ranges, 3, 5));
    assert!(order.takes(&ranges, 3, 5));
    assert!(!order.takes(&ranges, 3, 4));
    assert!(order.takes(&ranges, 3, 27));
    assert!(!ScanOrder::default().takes(&ranges, 3, 28));

    let ranges = [(0, 18), (21, 6), (9, 30)];
    let mut order = ScanOrder::default();
    assert!(order.takes(&ranges, 3, 15));
    assert!(order.takes(&ranges, 3, 9));
    assert!(!order.takes(&ranges, 3, 8));
    assert!(!order.takes(&ranges, 3, 5));
    assert!(!ScanOrder::default().takes(&[(0, 3), (21, 6)], 3, 21));

    let mut order = ScanOrder::default();
    assert!(order.takes(&[(5, 6), (10, 20)], 3, 5));
    assert!(!order.takes(&[(5, 6), (10, 20)], 3, 6));
    assert!(order.takes(&[(5, 6), (10, 20)], 3, 10));
}
