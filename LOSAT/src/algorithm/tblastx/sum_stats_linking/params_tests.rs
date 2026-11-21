use super::compute_avg_query_length_ncbi;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:200-215 (the contexts
// of the strand not searched have length 0) and c++/src/algo/blast/core/blast_parameters.c:
// 1023-1026 (the average over every context)
#[test]
fn average_query_length_counts_the_contexts_of_the_strand_not_searched() {
    // 300 bases: 100, 99 and 99 residues a strand. Both strands: offsets 0, 101, 201, 301, 402
    // and 502, (502 + 99 - 1) / 6 = 100.
    assert_eq!(compute_avg_query_length_ncbi(&[300], (true, true)), 100);
    // Plus: contexts 3 to 5 have no letters, at offset 301: (301 + 0 - 1) / 6 = 50.
    assert_eq!(compute_avg_query_length_ncbi(&[300], (true, false)), 50);
    // Minus: contexts 0 to 2 at offset 0, then 0, 101 and 201: (201 + 99 - 1) / 6 = 49.
    assert_eq!(compute_avg_query_length_ncbi(&[300], (false, true)), 49);
    // A second query of 30 bases: its empty contexts start after the first query's last.
    assert_eq!(compute_avg_query_length_ncbi(&[300, 30], (true, true)), 55);
    assert_eq!(compute_avg_query_length_ncbi(&[300, 30], (true, false)), 27);
    assert_eq!(compute_avg_query_length_ncbi(&[300, 30], (false, true)), 27);
    // No letters in any context: (0 + 0 - 1) / 6, 0 in C and in Rust.
    assert_eq!(compute_avg_query_length_ncbi(&[2], (false, true)), 0);
    // One query of 1,200 bases and 300 of 3 bases, one strand: 1,200 + 300 letters and the
    // separators in 1,806 contexts, an average of 0 (both strands: 1).
    let mut lengths = vec![1200];
    lengths.extend([3; 300]);
    assert_eq!(compute_avg_query_length_ncbi(&lengths, (true, false)), 0);
    assert_eq!(compute_avg_query_length_ncbi(&lengths, (false, true)), 0);
    assert_eq!(compute_avg_query_length_ncbi(&lengths, (true, true)), 1);
}
