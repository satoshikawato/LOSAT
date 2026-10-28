//! The hit lists of one query batch, filled inside its search: NCBI's HSP collector
//! (`hspfilter_collector.c`), which keeps at most `-max_target_seqs` subjects per query while the
//! subjects are searched. Before it, every HSP of every (query, subject) pair of the batch was
//! kept until the run's last report ordered them (`report::final_hit_order` over all batches), so
//! a `-db` search of many queries against a database that every query hits held the uncapped HSPs
//! of the whole run, about 190 bytes each (tblastx `-db 18S_fungal_sequences` x 50 of its records:
//! 3.7 GB, NCBI 0.43 GB; `blastdb-ids-20261021/R/REPORT_tbx_db_memory.md`).
//!
//! NCBI reference: c++/src/algo/blast/core/blast_engine.c:1543
//! ```c
//!             BlastHSPStreamWrite(hsp_stream, &hsp_list);
//! ```
//! NCBI reference: c++/src/algo/blast/core/hspfilter_collector.c:129-150
//! ```c
//!       for (index = 0; index < results->num_queries; index++) {
//!          if (hsp_list_array[index]) {
//!             if (!results->hitlist_array[index]) {
//!                results->hitlist_array[index] =
//!                   Blast_HitListNew(params->prelim_hitlist_size);
//!             }
//!             Blast_HitListUpdate(results->hitlist_array[index],
//!                                 hsp_list_array[index]);
//! ```
//! NCBI reference: c++/src/algo/blast/core/hspfilter_collector.c:337-338
//! ```c
//!    retval->prelim_hitlist_size = GetPrelimHitlistSize(hit_options->hitlist_size,
//!                                     compositionBasedStats, gapped_calculation);
//! ```
//! NCBI reference: c++/src/algo/blast/core/blast_hits.c:44-71 (`GetPrelimHitlistSize`): without
//! gapping or composition-based statistics the preliminary hit list has the final size, so
//! tblastx's lists are final when the batch's search ends; its traceback changes nothing.
//!
//! The subjects' HSP lists enter each query's hit list in OID order, as NCBI's single-threaded
//! subject loop writes them. With threads a subject's list comes when its search ends, so the
//! collector holds the lists of the subjects past the next OID until that one comes
//! (`HitCollector::push`); the loops hand out subjects in ascending order, so a few wait at a time.
//! Every subject of the loop sends one message, with or without HSPs, so the next OID always
//! comes. A query's lists are compared by best e-value, first HSP's score and OID
//! (`blastn::hsp::compare_hsp_lists`), a total order when no e-value is NaN, under which the kept
//! lists would not depend on the order they come in; the OID order makes the result that of
//! `report::final_hit_order` by construction (the same lists in the same order), NaN or not.
//!
//! With `-culling_limit` NCBI's culling writer replaces the collector (`setup_factory.cpp:330-341`);
//! the batch's HSPs are then kept until its end and culled query by query
//! (`report::culled_hit_order`), as before.

use std::collections::BTreeMap;

use crate::algorithm::blastn::hsp::{HitList, HitListEntry};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::core::blast_hits::trim_by_max_hsps;

use super::args::TblastxArgs;
use super::blast_engine::TblastxHsp;
use super::report::{self, HitFilters};

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
// ```c
// typedef struct BlastHSPList {
//    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
//    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
//                       Set to 0 if not applicable */
//    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
//    Int4 hspcnt; /**< Number of HSPs saved */
//    ...
//    double best_evalue; /**< Smallest e-value for HSPs in this list. Filled after
//                           e-values are calculated. Necessary because HSPs are
//                           sorted by score, but highest scoring HSP may not have
//                           the lowest e-value if sum statistics is used. */
// } BlastHSPList;
// ```
/// The final HSPs of one query and one subject, as the hit list holds them.
pub(super) struct TblastxHspList {
    pub(super) oid: u32,
    pub(super) hsps: Vec<TblastxHsp>,
    pub(super) best_evalue: f64,
}

impl HitListEntry for TblastxHspList {
    fn oid(&self) -> u32 {
        self.oid
    }

    fn hsp_count(&self) -> usize {
        self.hsps.len()
    }

    fn best_evalue(&self) -> f64 {
        self.best_evalue
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1739-1749
    // ```c
    // static double
    // s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
    // {
    //     int index = 0;
    //     double best_evalue = (double) INT4_MAX;
    //
    //     for (index=0; index<hsp_list->hspcnt; index++)
    //        best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
    //
    //     return best_evalue;
    // }
    // ```
    fn update_best_evalue(&mut self) {
        let mut best = i32::MAX as f64;
        for hsp in &self.hsps {
            if hsp.hit.e_value < best {
                best = hsp.hit.e_value;
            }
        }
        self.best_evalue = best;
    }

    fn first_score(&self) -> Option<i32> {
        self.hsps.first().map(|hsp| hsp.hit.raw_score)
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1452-1454
    // ```c
    //       qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //             s_EvalueCompareHSPs);
    // ```
    fn sort_by_evalue(&mut self) {
        self.hsps
            .sort_by(|a, b| crate::common::evalue_compare_hsps(&a.hit, &b.hit));
    }
}

/// The hit lists of the queries of one batch, by the query's index in the batch: filled one
/// subject's list at a time (`add`), emptied in report order (`finish`). The steps and their
/// NCBI references are `report::final_hit_order`'s, which runs on them.
pub(super) struct QueryLists {
    lists: Vec<Option<HitList<TblastxHspList>>>,
    hitlist_size: usize,
}

impl QueryLists {
    pub(super) fn new(hitlist_size: usize) -> Self {
        Self {
            lists: Vec::new(),
            hitlist_size,
        }
    }

    /// Subject `oid`'s HSPs of query `q_idx`, in the order of the subject's list (the state
    /// after linking): sorted by score, reaped by query coverage (`-qcov_hsp_perc`; a list left
    /// empty is not added), then `Blast_HitListUpdate` on the query's hit list.
    pub(super) fn add(
        &mut self,
        q_idx: u32,
        oid: u32,
        mut hsps: Vec<TblastxHsp>,
        filters: &HitFilters<'_>,
    ) {
        hsps.sort_by(|a, b| crate::common::score_compare_hsps(&a.hit, &b.hit));
        let context_lengths = filters.context_lengths(q_idx);
        if filters.reap_by_query_coverage(&mut hsps, &context_lengths) {
            return;
        }
        let index = q_idx as usize;
        if self.lists.len() <= index {
            self.lists.resize_with(index + 1, || None);
        }
        let hitlist_size = self.hitlist_size;
        self.lists[index]
            .get_or_insert_with(|| HitList::new(hitlist_size))
            .update(TblastxHspList {
                oid,
                hsps,
                best_evalue: 0.0,
            });
    }

    /// The HSPs of every query with a hit list, in report order: per query, each list cut to
    /// `-max_hsps`, the lists sorted by e-value, each list's HSPs sorted by e-value.
    pub(super) fn finish(self, filters: &HitFilters<'_>) -> Vec<TblastxHsp> {
        let mut ordered = Vec::new();
        for mut hit_list in self.lists.into_iter().flatten() {
            if filters.max_hsps > 0 {
                for list in &mut hit_list.hsplist_array {
                    trim_by_max_hsps(&mut list.hsps, filters.max_hsps);
                }
            }
            hit_list.sort_by_evalue();
            for mut list in hit_list.hsplist_array {
                HitListEntry::sort_by_evalue(&mut list);
                ordered.extend(list.hsps);
            }
        }
        ordered
    }
}

/// What the batch's search holds in place of NCBI's collector or culling writer.
enum Lists {
    /// The queries' hit lists (the default).
    Collector(QueryLists),
    /// `-culling_limit`: the batch's HSPs after the e-value reap, for `report::culled_hit_order`.
    Culling(Vec<TblastxHsp>),
}

/// Where the subject loop of a query batch's search puts each subject's final HSPs (`push`,
/// one message per subject of the loop), and what the batch's `hits` are made of (`finish`): the
/// queries' hit lists in report order, or, with `-culling_limit`, the culled lists.
pub(crate) struct HitCollector {
    lists: Lists,
    hitlist_size: usize,
    culling_limit: i32,
    max_hsps: usize,
    query_cov_hsp_perc: f64,
    /// The batch's queries' nucleotide lengths, read by the coverage reap and the culling
    /// (empty when neither applies).
    query_lengths: Vec<usize>,
    /// `-evalue`: an HSP above it is removed before the lists.
    evalue_threshold: f64,
    /// The OID the lists take next, and the subjects past it that have come (`push`).
    next: usize,
    waiting: BTreeMap<usize, Vec<TblastxHsp>>,
}

impl HitCollector {
    /// The collector of a batch of `queries` searched with `args`.
    pub(crate) fn new(args: &TblastxArgs, queries: &[FastaRecord]) -> Self {
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2960-2962,2978
        // ```c
        //     	if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
        //     		hitlist_size = args[kArgMaxTargetSequences].AsInteger();
        //     	}
        //     ...
        //     opt.SetHitlistSize(hitlist_size);
        // ```
        // The hit list size is `args.report_counts.hitlist_size` (`blastinput::report_counts`).
        // NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.hpp:1320-1324
        // ```c
        // CBlastOptionsLocal::SetCullingLimit(int s)
        // {
        //     if (s <= 0) {
        //         return;
        //     }
        // ```
        // A culling limit above 0 installs NCBI's culling writer and pipe
        // (`report::culled_hit_order`).
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:313-319
        // ```c
        //     if (args.Exist(kArgQueryCovHspPerc) && args[kArgQueryCovHspPerc]) {
        //         opt.SetQueryCovHspPerc(args[kArgQueryCovHspPerc].AsDouble());
        //     }
        //
        //     if (args.Exist(kArgMaxHSPsPerSubject) && args[kArgMaxHSPsPerSubject]) {
        //         opt.SetMaxHspsPerSubject(args[kArgMaxHSPsPerSubject].AsInteger());
        //     }
        // ```
        let culling_limit = args.culling_limit();
        let query_cov_hsp_perc = args.qcov_hsp_perc.unwrap_or(0.0);
        let query_lengths = if culling_limit > 0 || query_cov_hsp_perc > 0.0 {
            queries.iter().map(|query| query.seq().len()).collect()
        } else {
            Vec::new()
        };
        Self::with_options(
            args.report_counts.hitlist_size,
            culling_limit,
            args.max_hsps.unwrap_or(0),
            query_cov_hsp_perc,
            query_lengths,
            args.evalue,
        )
    }

    /// `new` with the options read: `-max_target_seqs`, `-culling_limit` (0: none), `-max_hsps`
    /// (0: none), `-qcov_hsp_perc` (0: none; `query_lengths` are then the queries' nucleotide
    /// lengths) and `-evalue`.
    pub(super) fn with_options(
        hitlist_size: usize,
        culling_limit: i32,
        max_hsps: usize,
        query_cov_hsp_perc: f64,
        query_lengths: Vec<usize>,
        evalue_threshold: f64,
    ) -> Self {
        let lists = if culling_limit > 0 {
            Lists::Culling(Vec::new())
        } else {
            Lists::Collector(QueryLists::new(hitlist_size))
        };
        Self {
            lists,
            hitlist_size,
            culling_limit,
            max_hsps,
            query_cov_hsp_perc,
            query_lengths,
            evalue_threshold,
            next: 0,
            waiting: BTreeMap::new(),
        }
    }

    /// Subject `s_idx`'s final HSPs (of every query of the batch, in the order of its list;
    /// none when it has none). The subjects reach the lists in OID order: one past the next
    /// OID waits here until that one has come.
    pub(crate) fn push(&mut self, s_idx: usize, hits: Vec<TblastxHsp>) {
        if s_idx != self.next {
            self.waiting.insert(s_idx, hits);
            return;
        }
        self.add_subject(s_idx, hits);
        self.next += 1;
        while let Some(hits) = self.waiting.remove(&self.next) {
            self.add_subject(self.next, hits);
            self.next += 1;
        }
    }

    // NCBI reference: c++/src/algo/blast/core/blast_hits.c:1988-1996
    // ```c
    //    cutoff = hit_options->expect_value;
    // ...
    //       if (hsp->evalue > cutoff) {
    // ```
    /// Subject `oid`'s HSPs into the lists: those above `-evalue` removed first (a NaN cutoff,
    /// `-evalue -nan`, keeps every HSP), then query by query, each query's HSPs in the order of
    /// the subject's list.
    fn add_subject(&mut self, oid: usize, mut hits: Vec<TblastxHsp>) {
        hits.retain(|hsp| !(hsp.hit.e_value > self.evalue_threshold));
        if hits.is_empty() {
            return;
        }
        let oid = u32::try_from(oid).expect("TBLASTX subject index out of range");
        match &mut self.lists {
            Lists::Culling(all) => all.extend(hits),
            Lists::Collector(lists) => {
                let filters = HitFilters {
                    max_hsps: self.max_hsps,
                    query_cov_hsp_perc: self.query_cov_hsp_perc,
                    query_lengths: &self.query_lengths,
                };
                hits.sort_by_key(|hsp| hsp.hit.q_idx);
                let mut rest = hits;
                while let Some(first) = rest.first() {
                    let q_idx = first.hit.q_idx;
                    let count = rest.iter().take_while(|hsp| hsp.hit.q_idx == q_idx).count();
                    let tail = rest.split_off(count);
                    lists.add(q_idx, oid, std::mem::replace(&mut rest, tail), &filters);
                }
            }
        }
    }

    /// The batch's HSPs in report order (`report::final_hit_order`'s, or `culled_hit_order`'s
    /// with `-culling_limit`). A subject still waiting (none when every subject of the loop
    /// sent) enters first, in OID order.
    pub(crate) fn finish(mut self) -> Vec<TblastxHsp> {
        for (oid, hits) in std::mem::take(&mut self.waiting) {
            self.add_subject(oid, hits);
        }
        let filters = HitFilters {
            max_hsps: self.max_hsps,
            query_cov_hsp_perc: self.query_cov_hsp_perc,
            query_lengths: &self.query_lengths,
        };
        match self.lists {
            Lists::Collector(lists) => lists.finish(&filters),
            Lists::Culling(hits) => report::culled_hit_order(
                hits,
                self.hitlist_size,
                self.culling_limit,
                &self.query_lengths,
                &filters,
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Hit;

    fn hsp(q_idx: u32, s_idx: u32, raw_score: i32, e_value: f64, offset: usize) -> TblastxHsp {
        TblastxHsp {
            hit: Hit {
                identity: 0.0,
                length: 10,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 30,
                s_start: 1,
                s_end: 30,
                e_value,
                bit_score: raw_score as f64 / 2.0,
                num_ident: 0,
                query_frame: 1,
                query_length: 0,
                q_idx,
                s_idx,
                raw_score,
                sort_query_offset: offset,
                sort_query_end: offset + 10,
                sort_subject_offset: offset,
                sort_subject_end: offset + 10,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 0,
            },
            subject_frame: 1,
            num: 1,
        }
    }

    fn key(hsps: &[TblastxHsp]) -> Vec<(u32, u32, i32, f64, usize)> {
        hsps.iter()
            .map(|hsp| {
                (
                    hsp.hit.q_idx,
                    hsp.hit.s_idx,
                    hsp.hit.raw_score,
                    hsp.hit.e_value,
                    hsp.hit.sort_query_offset,
                )
            })
            .collect()
    }

    // The final HSPs of `subjects` subjects for two queries, as the subject loop sends them (one
    // list per subject, the queries' HSPs interleaved, in score order after linking). Every
    // seventh subject has none; most best e-values tie (0.0, as 18S rRNA against its database),
    // and in many lists the HSP with the best e-value is not the one with the highest score, so
    // the comparison of the lists reads the first HSP differently before and after the hit list
    // fills. A few HSPs are above the e-value threshold 1e-5.
    fn subject_lists(subjects: usize) -> Vec<(usize, Vec<TblastxHsp>)> {
        (0..subjects)
            .map(|s| {
                let oid = s as u32;
                if s % 7 == 6 {
                    return (s, Vec::new());
                }
                let best = if s % 3 == 0 {
                    0.0
                } else {
                    1e-6 * (1 + s % 5) as f64
                };
                let mut hsps = vec![
                    hsp(1, oid, 130 + (s % 4) as i32, 1e-3, 40),
                    hsp(0, oid, 120 + (s % 4) as i32, 1e-3, 0),
                    hsp(0, oid, 100 + (s % 6) as i32, best, 20),
                    hsp(1, oid, 110 + (s % 6) as i32, best * 0.5, 60),
                ];
                if s % 11 == 10 {
                    hsps.push(hsp(0, oid, 50, 1e-2, 80));
                }
                if s % 2 == 1 {
                    hsps.push(hsp(0, oid, 100 + (s % 6) as i32, best, 10));
                }
                (s, hsps)
            })
            .collect()
    }

    fn reference(
        lists: &[(usize, Vec<TblastxHsp>)],
        hitlist_size: usize,
        filters: &HitFilters<'_>,
    ) -> Vec<TblastxHsp> {
        let hits: Vec<TblastxHsp> = lists
            .iter()
            .flat_map(|(_, hsps)| hsps.iter().cloned())
            .filter(|hsp| !(hsp.hit.e_value > 1e-5))
            .collect();
        report::final_hit_order(hits, hitlist_size, filters)
    }

    // A permutation of 0..n (7919 is a prime above every n used here).
    fn shuffled(n: usize) -> Vec<usize> {
        (0..n).map(|i| i * 7919 % n).collect()
    }

    // NCBI c++/src/algo/blast/core/hspfilter_collector.c:129-150 and blast_hits.c:3243-3297
    // (Blast_HitListUpdate): the subjects enter the hit lists in OID order whatever order their
    // searches end in, and the lists are those of the whole batch's HSPs ordered at once
    // (`final_hit_order`), below, at and above the hit list size.
    #[test]
    fn subjects_in_any_order_give_final_hit_orders_lists() {
        for subjects in [12, 499, 500, 501, 520] {
            let lists = subject_lists(subjects);
            let want = reference(&lists, 500, &HitFilters::default());
            for order in [(0..subjects).collect::<Vec<_>>(), shuffled(subjects)] {
                let mut collector = HitCollector::with_options(500, 0, 0, 0.0, Vec::new(), 1e-5);
                for &s in &order {
                    collector.push(s, lists[s].1.clone());
                }
                let got = collector.finish();
                assert_eq!(key(&got), key(&want), "{subjects} subjects");
            }
        }
    }

    // `-max_hsps` cuts each kept list, `-qcov_hsp_perc` reaps a subject's HSPs before its list
    // enters the hit list: as `final_hit_order` applies them.
    #[test]
    fn filters_apply_as_in_final_hit_order() {
        let subjects = 503;
        let lists = subject_lists(subjects);
        let query_lengths = vec![300usize, 90];
        let filters = HitFilters {
            max_hsps: 2,
            query_cov_hsp_perc: 25.0,
            query_lengths: &query_lengths,
        };
        let want = reference(&lists, 500, &filters);
        assert!(!want.is_empty());
        let mut collector =
            HitCollector::with_options(500, 0, 2, 25.0, query_lengths.clone(), 1e-5);
        for s in shuffled(subjects) {
            collector.push(s, lists[s].1.clone());
        }
        assert_eq!(key(&collector.finish()), key(&want));
    }

    // `-culling_limit`: the batch's HSPs are culled at its end (`culled_hit_order`).
    #[test]
    fn culling_limit_culls_the_batch_at_its_end() {
        let subjects = 60;
        let lists = subject_lists(subjects);
        let query_lengths = vec![300usize, 90];
        let hits: Vec<TblastxHsp> = lists
            .iter()
            .flat_map(|(_, hsps)| hsps.iter().cloned())
            .filter(|hsp| !(hsp.hit.e_value > 1e-5))
            .collect();
        let want = report::culled_hit_order(
            hits,
            500,
            1,
            &query_lengths,
            &HitFilters {
                max_hsps: 0,
                query_cov_hsp_perc: 0.0,
                query_lengths: &query_lengths,
            },
        );
        assert!(!want.is_empty());
        let mut collector = HitCollector::with_options(500, 1, 0, 0.0, query_lengths.clone(), 1e-5);
        for s in shuffled(subjects) {
            collector.push(s, lists[s].1.clone());
        }
        assert_eq!(key(&collector.finish()), key(&want));
    }

    // A subject that never sends (no subject of the loop today) does not lose the ones past it:
    // they enter at the end, in OID order.
    #[test]
    fn subjects_past_a_missing_one_enter_at_the_end() {
        let lists: Vec<(usize, Vec<TblastxHsp>)> = subject_lists(30)
            .into_iter()
            .map(|(s, hsps)| {
                if s == 4 || s == 17 {
                    (s, Vec::new())
                } else {
                    (s, hsps)
                }
            })
            .collect();
        let want = reference(&lists, 500, &HitFilters::default());
        let mut collector = HitCollector::with_options(500, 0, 0, 0.0, Vec::new(), 1e-5);
        for s in shuffled(30) {
            if s != 4 && s != 17 {
                collector.push(s, lists[s].1.clone());
            }
        }
        assert_eq!(key(&collector.finish()), key(&want));
    }
}
