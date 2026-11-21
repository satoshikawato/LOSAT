use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1497-1499 (a third of the query's letters), :1381-1384 and
// showalign.cpp:3162-3185 with tabular.cpp:737-747 (ungapped: the query's letters, or the
// span's nucleotides for a subject's single HSP). Values of tab.bx.cov.6 (qn2: 135 letters).
#[test]
fn hsp_coverage_follows_ncbis_seq_aligns() {
    assert_eq!(query_coverage_per_hsp(true, 2, 6, 18, 135, 135), 13);
    assert_eq!(query_coverage_per_hsp(true, 1, 15, 45, 135, 135), 33);
    assert_eq!(query_coverage_per_hsp(false, 2, 6, 18, 135, 135), 4);
    assert_eq!(query_coverage_per_hsp(false, 1, 6, 18, 135, 135), 13);
    assert_eq!(query_coverage_per_hsp(true, 1, 140, 420, 420, 420), 100);
    // With `-query_loc`: the interval for the engine's coverage, the range as given for the
    // formatter's.
    assert_eq!(query_coverage_per_hsp(true, 1, 15, 45, 135, 300), 33);
    assert_eq!(query_coverage_per_hsp(false, 2, 6, 18, 135, 300), 4);
    assert_eq!(query_coverage_per_hsp(false, 1, 6, 18, 135, 300), 6);
}
