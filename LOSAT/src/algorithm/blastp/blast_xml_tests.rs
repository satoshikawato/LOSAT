//! Unit tests of BLASTP's BLAST XML adapter: PHI-BLAST's results (`XmlResult`).

use super::*;

// Line PS, PS15b-4: the results of one query's pattern occurrences join the held report of the
// first (a `-subject` search writes them as one report, NCBI BLAST+ 2.17.0 `PS/ps15b4/ncbi_phi`
// `h1.tk2.sub.16`: two reports of three `Search`es); another query's starts its own. Without the
// run's scopes no result shows another's label.
#[test]
fn a_phi_result_joins_the_held_report_of_its_query() {
    let results = [
        XmlResult {
            query: 0,
            warned: 1,
        },
        XmlResult {
            query: 0,
            warned: 1,
        },
        XmlResult {
            query: 1,
            warned: 2,
        },
    ];
    let shows: Vec<_> = (0..results.len())
        .map(|index| result_show(None, &results, index))
        .collect();
    assert_eq!(
        shows.iter().map(|show| show.joins_held).collect::<Vec<_>>(),
        [false, true, false]
    );
    assert!(shows.iter().all(|show| show.label.is_none()));
}
