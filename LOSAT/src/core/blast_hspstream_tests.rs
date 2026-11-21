use super::{cbs_trim_count, cbs_trim_sorts};

/// The count a hit list of the best E-values `evalues` (in E-value order) keeps after the
/// trim, as a caller runs it: the test, the sort (here the order given) and the cut.
fn trimmed(evalues: &[f64], hitlist_size: usize) -> usize {
    if !cbs_trim_sorts(evalues.len(), hitlist_size) {
        return evalues.len();
    }
    cbs_trim_count(evalues.len(), evalues.len(), hitlist_size, |i| evalues[i])
        .unwrap_or(evalues.len())
}

/// `count` E-values: `low` for the first `split`, then `high`.
fn list(count: usize, split: usize, low: f64, high: f64) -> Vec<f64> {
    (0..count)
        .map(|i| if i < split { low } else { high })
        .collect()
}

// No NCBI unit test covers BlastHSPCBSStreamClose; the cases follow its arithmetic
// (blast_hspstream.c:829-856) at the edges of each branch.
#[test]
fn the_list_is_cut_at_the_first_hundred_past_the_limit() {
    // 1,200 lists are not more than MIN_BUF_SZ + 100 = 1,200: nothing happens.
    assert!(!cbs_trim_sorts(1200, 500));
    assert_eq!(trimmed(&list(1200, 700, 1e-30, 1e-5), 500), 1200);
    // 1,201: E[500] = 1e-30, mag -30 -> max(-27, -20) = -20, limit 9.9e-20; E[1100] =
    // 1e-5 passes it: the list keeps 1,100.
    assert!(cbs_trim_sorts(1201, 500));
    assert_eq!(trimmed(&list(1201, 700, 1e-30, 1e-5), 500), 1100);
    // The first index past the limit is 1,250: the cut is at 1,300 (steps of 100), and
    // 1,300 is max_index of a 1,301-list, which the loop does not reach.
    assert_eq!(trimmed(&list(1400, 1250, 1e-30, 1e-5), 500), 1300);
    assert_eq!(trimmed(&list(1301, 1250, 1e-30, 1e-5), 500), 1301);
    // A larger hit list size moves both indices: REF 700, MIN_BUF_SZ 1,300.
    assert!(!cbs_trim_sorts(1400, 700));
    assert_eq!(trimmed(&list(1500, 900, 1e-30, 1e-5), 700), 1300);
    // A hit list size below 500 counts as 500.
    assert_eq!(trimmed(&list(1201, 700, 1e-30, 1e-5), 10), 1100);
}

#[test]
fn the_limit_follows_ncbi_s_integer_magnitudes() {
    // E[500] = 1e-10: mag -10 -> -10/2 = -5, limit 9.9e-5; 1e-4 passes it, 9e-5 does not.
    assert_eq!(trimmed(&list(1300, 1000, 1e-10, 1e-4), 500), 1100);
    assert_eq!(trimmed(&list(1300, 1000, 1e-10, 9e-5), 500), 1300);
    // E[500] = 3e-25: log10 = -24.5 -> -24 (toward zero), -24*90/100 = -21 (toward zero),
    // max(-21, -14) = -14, limit 9.9e-14.
    assert_eq!(trimmed(&list(1300, 1000, 3e-25, 1e-13), 500), 1100);
    assert_eq!(trimmed(&list(1300, 1000, 3e-25, 9e-14), 500), 1300);
    // E[500] = 5: mag 0 -> 0, limit 3 * 5 = 15.
    assert_eq!(trimmed(&list(1300, 1000, 5.0, 16.0), 500), 1100);
    // E[500] = 0: mag stays -180 -> max(-162, -170) = -162; an E-value of 0 never cuts.
    assert_eq!(trimmed(&list(1300, 1300, 0.0, 0.0), 500), 1300);
    assert_eq!(trimmed(&list(1300, 1000, 0.0, 1e-160), 500), 1100);
}

#[test]
fn a_list_shortened_by_the_purge_is_read_to_its_sorted_end() {
    // 1,300 lists of which 150 were empty: the loop reads indices 1,100 and below 1,150.
    let evalues = list(1150, 1000, 1e-30, 1e-5);
    assert_eq!(cbs_trim_count(1300, 1150, 500, |i| evalues[i]), Some(1100));
    let late = list(1150, 1120, 1e-30, 1e-5);
    assert_eq!(cbs_trim_count(1300, 1150, 500, |i| late[i]), None);
}
