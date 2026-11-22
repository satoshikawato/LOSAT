//! The hit lists of one query batch, filled inside its search: NCBI's HSP collector
//! (`hspfilter_collector.c`), which keeps at most `-max_target_seqs` subjects per query while the
//! subjects are searched. Before it, every HSP of every (query, subject) pair of the batch was
//! kept until the run's last report ordered them (`report::final_hit_order` over all batches), so
//! a `-db` search of many queries against a database that every query hits held the uncapped HSPs
//! of the whole run, about 190 bytes each (tblastx `-db 18S_fungal_sequences` x 50 of its records:
//! 3.7 GB, NCBI 0.43 GB; `blastdb-ids-20261021/R/REPORT_tbx_db_memory.md`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1554
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:143-150
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:336-337
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:44-71 (`GetPrelimHitlistSize`): without
//! gapping or composition-based statistics the preliminary hit list has the final size, so
//! tblastx's lists are final when the batch's search ends; its traceback changes nothing.
//!
//! The subjects' HSP lists enter each query's hit list in OID order, as NCBI's single-threaded
//! subject loop writes them. With threads a subject's list comes when its search ends, so the
//! collector holds the lists of the subjects past the next OID until that one comes
//! (`HitCollector::push`). The default loop (LOSAT_X_TBXNEXT, `x_subject_next.rs`) hands out
//! subjects in ascending order, so a few wait at a time, more while one long subject is searched;
//! rayon's pieces (the loop without the switch) start far apart, so most of a batch waits until
//! its first piece ends, as the whole batch did before. Every subject of the loop sends one
//! message, with or without HSPs, so the next OID always comes.
//! A query's lists are compared by best e-value, first HSP's score and OID
//! (`blastn::hsp::compare_hsp_lists`), a total order when no e-value is NaN, under which the kept
//! lists would not depend on the order they come in; the OID order makes the result that of
//! `report::final_hit_order` by construction (the same lists in the same order), NaN or not.
//!
//! With `-culling_limit` NCBI's culling writer replaces the collector (`setup_factory.cpp:330-341`);
//! the batch's HSPs are then kept until its end and culled query by query
//! (`report::culled_hit_order`), as before. With `-best_hit_overhang` or `-best_hit_score_edge`
//! NCBI's best hit writer replaces it (`besthit::BestHitWriter`), filled in the same OID order.

use std::collections::BTreeMap;

use crate::algorithm::blastn::hsp::{HitList, HitListEntry};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::core::blast_hits::trim_by_max_hsps;

use crate::core::hspfilter_besthit::BestHitParams;

use super::args::TblastxArgs;
use super::besthit::BestHitWriter;
use super::blast_engine::TblastxHsp;
use super::report::{self, HitFilters};

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1739-1749
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1452-1454
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

    /// Subject `oid`'s HSPs of query `q_idx` (`subject_list`): `Blast_HitListUpdate` on the
    /// query's hit list.
    pub(super) fn add(&mut self, q_idx: u32, oid: u32, hsps: Vec<TblastxHsp>) {
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
    /// `-max_hsps`, the lists sorted by e-value (`sorts_results`) or in the traceback's order
    /// where NCBI does not sort them (`HitList::put_in_traceback_order`: a database length of 0 or
    /// of 2^63 or more; the preliminary list's heap order is not the traceback's).
    pub(super) fn finish(self, filters: &HitFilters<'_>, sorts_results: bool) -> Vec<TblastxHsp> {
        let mut ordered = Vec::new();
        let hitlist_size = self.hitlist_size;
        for mut hit_list in self.lists.into_iter().flatten() {
            if !sorts_results {
                hit_list.put_in_traceback_order(hitlist_size);
            }
            if filters.max_hsps > 0 || filters.subject_besthit {
                for list in &mut hit_list.hsplist_array {
                    trim_by_max_hsps(&mut list.hsps, filters.max_hsps);
                    filters.subject_besthit(&mut list.hsps);
                }
            }
            if sorts_results {
                hit_list.sort_by_evalue();
            }
            for mut list in hit_list.hsplist_array {
                HitListEntry::sort_by_evalue(&mut list);
                ordered.extend(list.hsps);
            }
        }
        ordered
    }
}

/// Subject `oid`'s HSPs of query `q_idx` as NCBI's preliminary stage writes them, in the order
/// of the subject's list (the state after linking): sorted by score and reaped by query coverage
/// (`-qcov_hsp_perc`); `None` when none is left (NCBI writes no empty list).
pub(super) fn subject_list(
    q_idx: u32,
    mut hsps: Vec<TblastxHsp>,
    filters: &HitFilters<'_>,
) -> Option<Vec<TblastxHsp>> {
    hsps.sort_by(|a, b| crate::common::score_compare_hsps(&a.hit, &b.hit));
    let context_lengths = filters.context_lengths(q_idx);
    if filters.reap_by_query_coverage(&mut hsps, &context_lengths) {
        return None;
    }
    Some(hsps)
}

/// What the batch's search holds in place of NCBI's collector or culling writer.
enum Lists {
    /// The queries' hit lists (the default).
    Collector(QueryLists),
    /// `-best_hit_overhang` / `-best_hit_score_edge`: NCBI's best hit writer, whose `finish`
    /// runs the pipe after the traceback (`besthit`).
    BestHit(BestHitWriter),
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
    /// `-subject_besthit` (`HitFilters::subject_besthit`).
    subject_besthit: bool,
    /// The batch's queries' nucleotide lengths, read by the coverage reap, the culling and the
    /// best hit filter (empty when none applies).
    query_lengths: Vec<usize>,
    /// `-evalue`: an HSP above it is removed before the lists.
    evalue_threshold: f64,
    /// Whether NCBI's traceback sorts the results at its end (a positive `BlastSeqSrcGetTotLen`:
    /// `SubjectSetNumbers::link_total`).
    sorts_results: bool,
    /// The OID the lists take next, and the subjects past it that have come (`push`).
    next: usize,
    waiting: BTreeMap<usize, Vec<TblastxHsp>>,
}

impl HitCollector {
    /// The collector of a batch of `queries` searched with `args`.
    pub(crate) fn new(args: &TblastxArgs, queries: &[FastaRecord], sorts_results: bool) -> Self {
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2960-2962,2978
        // The hit list size is `args.report_counts.hitlist_size` (`blastinput::report_counts`).
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1320-1324
        // A culling limit above 0 installs NCBI's culling writer and pipe
        // (`report::culled_hit_order`).
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:313-319
        let culling_limit = args.culling_limit();
        let query_cov_hsp_perc = args.qcov_hsp_perc.unwrap_or(0.0);
        let best_hit = args.best_hit();
        let query_lengths = if culling_limit > 0 || query_cov_hsp_perc > 0.0 || best_hit.is_some() {
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
        .sorting_results(sorts_results)
        .with_subject_besthit(args.subject_besthit)
        .with_best_hit(best_hit)
    }

    /// The collector with NCBI's best hit writer in place of its HSP collector when `best_hit`
    /// is given (`-best_hit_overhang`, `-best_hit_score_edge`; never with `-culling_limit`,
    /// which the parser refuses). `query_lengths` must then be the queries' nucleotide lengths.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329
    pub(super) fn with_best_hit(self, best_hit: Option<BestHitParams>) -> Self {
        match best_hit {
            Some(params) if self.culling_limit <= 0 => Self {
                lists: Lists::BestHit(BestHitWriter::new(
                    params,
                    self.hitlist_size,
                    &self.query_lengths,
                )),
                ..self
            },
            _ => self,
        }
    }

    /// The collector with `-subject_besthit` (`HitCollector::subject_besthit`).
    pub(super) fn with_subject_besthit(self, subject_besthit: bool) -> Self {
        Self {
            subject_besthit,
            ..self
        }
    }

    /// The collector with `sorts_results` (`HitCollector::sorts_results`).
    pub(super) fn sorting_results(self, sorts_results: bool) -> Self {
        Self {
            sorts_results,
            ..self
        }
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
            subject_besthit: false,
            query_lengths,
            evalue_threshold,
            sorts_results: true,
            next: 0,
            waiting: BTreeMap::new(),
        }
    }

    /// Subject `s_idx`'s final HSPs (of every query of the batch, in the order of its list;
    /// none when it has none). The subjects reach the lists in OID order: one past the next
    /// OID waits here until that one has come.
    pub(crate) fn push(&mut self, s_idx: usize, hits: Vec<TblastxHsp>) {
        debug_assert!(
            s_idx >= self.next && !self.waiting.contains_key(&s_idx),
            "subject {s_idx} sent twice"
        );
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1988-1996
    /// Subject `oid`'s HSPs into the lists: those above `-evalue` removed first (a NaN cutoff,
    /// `-evalue -nan`, keeps every HSP), then query by query, each query's HSPs in the order of
    /// the subject's list.
    fn add_subject(&mut self, oid: usize, mut hits: Vec<TblastxHsp>) {
        #[expect(
            clippy::neg_cmp_op_on_partial_ord,
            reason = "NCBI removes an HSP only when evalue > cutoff (blast_hits.c:1996); `<=` would drop NaN values"
        )]
        hits.retain(|hsp| !(hsp.hit.e_value > self.evalue_threshold));
        if hits.is_empty() {
            return;
        }
        let oid = u32::try_from(oid).expect("TBLASTX subject index out of range");
        match &mut self.lists {
            Lists::Culling(all) => all.extend(hits),
            lists => {
                let filters = HitFilters {
                    max_hsps: self.max_hsps,
                    query_cov_hsp_perc: self.query_cov_hsp_perc,
                    query_lengths: &self.query_lengths,
                    subject_besthit: self.subject_besthit,
                };
                // Each query's HSPs moved once into a list grown by `push` (a list is kept
                // until the batch ends; `split_off` would keep every list at the capacity of
                // the whole subject's).
                hits.sort_by_key(|hsp| hsp.hit.q_idx);
                let mut hits = hits.into_iter().peekable();
                while let Some(q_idx) = hits.peek().map(|hsp| hsp.hit.q_idx) {
                    let mut list = Vec::new();
                    while let Some(hsp) = hits.next_if(|hsp| hsp.hit.q_idx == q_idx) {
                        list.push(hsp);
                    }
                    let Some(list) = subject_list(q_idx, list, &filters) else {
                        continue;
                    };
                    match lists {
                        Lists::Collector(lists) => lists.add(q_idx, oid, list),
                        Lists::BestHit(writer) => writer.add(q_idx, oid, list),
                        Lists::Culling(_) => unreachable!("the culling lists take every HSP above"),
                    }
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
            subject_besthit: self.subject_besthit,
        };
        match self.lists {
            Lists::Collector(lists) => lists.finish(&filters, self.sorts_results),
            Lists::BestHit(writer) => writer.finish(&filters, self.sorts_results),
            Lists::Culling(hits) => report::culled_hit_order(
                hits,
                self.hitlist_size,
                self.culling_limit,
                &self.query_lengths,
                &filters,
                self.sorts_results,
            ),
        }
    }
}

#[cfg(test)]
#[path = "hit_collector_tests.rs"]
mod tests;
