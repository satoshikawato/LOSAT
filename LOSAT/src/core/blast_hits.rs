//! HSP Storage and Management
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c
//!
//! This module provides HSP storage, sorting, and management
//! functionality.
//!
//! During migration, this re-exports from the original locations. The filters of
//! `s_FilterBlastResults` that every program shares (`-max_hsps`, `-qcov_hsp_perc`) are here.

pub use crate::algorithm::tblastx::chaining::*;

/// The query coverage of an HSP in percent: its query extent (`offset` to `end`, end
/// exclusive, in the letters of its context) over the length of its context, plus 0.5 below
/// 99 (`Blast_HSPGetQueryCoverage`).
///
/// NCBI reference: c++/src/algo/blast/core/blast_hits.c:1034-1044
/// ```c
/// double Blast_HSPGetQueryCoverage(const BlastHSP* hsp, Int4 query_length)
/// {
/// 	double pct = 0;
///     if(query_length > 0) {
///             pct = 100.0 * (double) (hsp->query.end - hsp->query.offset)/ (double) query_length;
///             if(pct < 99)
///             	pct +=0.5;
///     }
///     return pct;
/// }
/// ```
pub fn hsp_query_coverage(offset: i64, end: i64, query_length: i64) -> f64 {
    let mut pct = 0.0;
    if query_length > 0 {
        pct = 100.0 * (end - offset) as f64 / query_length as f64;
        if pct < 99.0 {
            pct += 0.5;
        }
    }
    pct
}

/// Whether an HSP covers less than `min_pct` percent of its query context
/// (`Blast_HSPQueryCoverageTest`: the HSP is deleted).
///
/// NCBI reference: c++/src/algo/blast/core/blast_hits.c:1046-1051
/// ```c
/// Boolean Blast_HSPQueryCoverageTest(BlastHSP* hsp,
///         						   double min_query_coverage_pct,
///         						   Int4 query_length)
/// {
///      double hsp_coverage = Blast_HSPGetQueryCoverage( hsp, query_length);
///      return (hsp_coverage < min_query_coverage_pct);
/// }
/// ```
pub fn hsp_query_coverage_test(min_pct: f64, offset: i64, end: i64, query_length: i64) -> bool {
    hsp_query_coverage(offset, end, query_length) < min_pct
}

/// `Blast_HSPListReapByQueryCoverage` over the HSPs of one list, in their order: the HSPs
/// whose query coverage (`hsp_query_coverage_test` of the `extent` of each: query offset,
/// query end and the length of its context) is below `min_pct` are removed. A `min_pct` of
/// 0 is NCBI's unset option and removes nothing. Returns whether an HSP was removed (the
/// caller then sets the list's best e-value again, and frees a list left empty).
///
/// NCBI reference: c++/src/algo/blast/core/blast_hits.c:2010-2045
/// ```c
///    if ((hsp_list == NULL) || (hsp_list->hspcnt == 0) || (hit_options->query_cov_hsp_perc == 0))
///       return 0;
///
///    hsp_array = hsp_list->hsp_array;
///    for (index = 0; index < hsp_list->hspcnt; index++) {
///       hsp = hsp_array[index];
///       ASSERT(hsp != NULL);
///       if ( Blast_HSPQueryCoverageTest(hsp, hit_options->query_cov_hsp_perc,
/// 			                          query_info->contexts[hsp->context].query_length)) {
///          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
///          removed = TRUE;
///       } else {
///          if (index > hsp_cnt)
///             hsp_array[hsp_cnt] = hsp_array[index];
///          hsp_cnt++;
///       }
///    }
///
///    hsp_list->hspcnt = hsp_cnt;
///
///    if (removed) {
///        hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
///    }
/// ```
pub fn reap_by_query_coverage<T>(
    hsps: &mut Vec<T>,
    min_pct: f64,
    extent: impl Fn(&T) -> (i64, i64, i64),
) -> bool {
    if hsps.is_empty() || min_pct == 0.0 {
        return false;
    }
    let count = hsps.len();
    hsps.retain(|hsp| {
        let (offset, end, query_length) = extent(hsp);
        !hsp_query_coverage_test(min_pct, offset, end, query_length)
    });
    hsps.len() != count
}

/// `Blast_TrimHSPListByMaxHsps`: the first `max_hsps` HSPs of a list, in its order, are
/// kept (0 keeps all). The list's best e-value is not set again.
///
/// NCBI reference: c++/src/algo/blast/core/blast_hits.c:2049-2069
/// ```c
/// Int2 Blast_TrimHSPListByMaxHsps(BlastHSPList* hsp_list,
///                                 const BlastHitSavingOptions* hit_options)
/// {
///    ...
///    if ((hsp_list == NULL) ||
/// 	   (hit_options->max_hsps_per_subject == 0) ||
/// 	   (hsp_list->hspcnt <= hit_options->max_hsps_per_subject))
///       return 0;
///
///    hsp_max = hit_options->max_hsps_per_subject;
///    hsp_array = hsp_list->hsp_array;
///    for (index = hsp_max; index < hsp_list->hspcnt; index++) {
///       hsp_array[index] = Blast_HSPFree(hsp_array[index]);
///    }
///
///    hsp_list->hspcnt = hsp_max;
///    return 0;
/// }
/// ```
pub fn trim_by_max_hsps<T>(hsps: &mut Vec<T>, max_hsps: usize) {
    if max_hsps == 0 || hsps.len() <= max_hsps {
        return;
    }
    hsps.truncate(max_hsps);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn query_coverage_adds_a_half_below_99() {
        // 95 of 1000 letters: 9.5 % + 0.5.
        assert_eq!(hsp_query_coverage(5, 100, 1000), 10.0);
        // 988 of 1000: 98.8 % + 0.5 = 99.3; 992 of 1000 stays 99.2.
        assert_eq!(hsp_query_coverage(0, 988, 1000), 98.8 + 0.5);
        assert_eq!(hsp_query_coverage(0, 992, 1000), 99.2);
        assert_eq!(hsp_query_coverage(0, 10, 0), 0.0);
    }

    #[test]
    fn reap_keeps_the_order_and_reports_removals() {
        let mut hsps = vec![(0, 50), (0, 95), (10, 100), (0, 9)];
        let extent = |&(offset, end): &(i64, i64)| (offset, end, 1000);
        assert!(!reap_by_query_coverage(&mut hsps, 0.0, extent));
        assert_eq!(hsps.len(), 4);
        // Coverage 5.5, 10.0, 9.5, 1.4: a limit of 9.5 keeps the second and the third.
        assert!(reap_by_query_coverage(&mut hsps, 9.5, extent));
        assert_eq!(hsps, vec![(0, 95), (10, 100)]);
        assert!(!reap_by_query_coverage(&mut hsps, 9.5, extent));
    }

    #[test]
    fn trim_keeps_the_first_hsps() {
        let mut hsps = vec![3, 1, 2];
        trim_by_max_hsps(&mut hsps, 0);
        assert_eq!(hsps, vec![3, 1, 2]);
        trim_by_max_hsps(&mut hsps, 3);
        assert_eq!(hsps, vec![3, 1, 2]);
        trim_by_max_hsps(&mut hsps, 2);
        assert_eq!(hsps, vec![3, 1]);
    }
}
