//! NCBI BLAST best hit filtering (`-best_hit_overhang`, `-best_hit_score_edge`): the lists of
//! `c++/src/algo/blast/core/hspfilter_besthit.c`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c
//!
//! Either option turns the filter on, the other taking its default (0.1). NCBI filters twice:
//! a writer in the preliminary stage in place of the HSP collector, whose overhang is 0, and a
//! pipe after the traceback with both values (`CSetupFactory::CreateHspWriter`,
//! `CreateHspPipe`). Each query keeps one list of HSPs ordered by their start on the plus
//! strand of the query; a new HSP is dropped when an HSP of the list contains it with a lower
//! e-value and a higher score per letter (by the score edge), and it removes the HSPs of the
//! list that it contains (by the overhang) with a higher e-value and a lower score per letter.
//!
//! This module holds no program type: the list never touches the payload `H` of an HSP. The
//! keys it reads (the range on the plus strand of the query, the score and the e-value:
//! `hsp->query`, `hsp->score`, `hsp->evalue` of `s_BlastHSPBestHitRun`) are made by the caller,
//! which also decides the query and how a minus context folds its range.
//!
//! NCBI's linked lists are vectors here, in the same order: a list is always ordered by its
//! `begin` (every insertion goes before the first entry whose `begin` is not smaller), so the
//! searches by `begin` are binary searches that find the node NCBI's walk finds.
//!
//! The overhangs are NCBI's `(Int4)` conversions of a double (`ncbi_int4_from_double`: an
//! overhang near 0.5 makes `2 * len * overhang / (1 - 2 * overhang)` larger than an `Int4`, which
//! the oracle's conversion makes `INT_MIN`, as for the X-drops, Owner 2026-11-17). The ranges and
//! the counts are `i64`: NCBI's `begin - overhang` with that `INT_MIN`, and
//! `prelim_hitlist_size * 2` past an `Int4`, are signed overflows (undefined in C, not reproduced,
//! NCBI parity DECISIONS 17). With `INT_MIN` the range of the HSPs to remove is empty, as it is in
//! the oracle, whose wrapped range lies below every HSP.

use std::cmp::Ordering;
use std::collections::HashMap;

use crate::core::blast_util::ncbi_int4_from_double;

/// NCBI reference (598d8ae6): c++/include/algo/blast/core/hspfilter_besthit.h:80-91
/// The defaults of the two options: an option that is not given takes this value when the
/// other is given (`CBlastOptionsLocal::SetBestHitOverhang`, `SetBestHitScoreEdge`).
pub(crate) const BEST_HIT_OVERHANG_DEFAULT: f64 = 0.1;
pub(crate) const BEST_HIT_SCORE_EDGE_DEFAULT: f64 = 0.1;

/// The values of the best hit filter of a search: present when `-best_hit_overhang` or
/// `-best_hit_score_edge` is given.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1354-1393
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3341-3346
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct BestHitParams {
    pub(crate) overhang: f64,
    pub(crate) score_edge: f64,
}

impl BestHitParams {
    /// The filter of the two options as given; `None` when neither is.
    pub(crate) fn from_options(overhang: Option<f64>, score_edge: Option<f64>) -> Option<Self> {
        if overhang.is_none() && score_edge.is_none() {
            return None;
        }
        Some(Self {
            overhang: overhang.unwrap_or(BEST_HIT_OVERHANG_DEFAULT),
            score_edge: score_edge.unwrap_or(BEST_HIT_SCORE_EDGE_DEFAULT),
        })
    }
}

/// What `s_BlastHSPBestHitRun` reads of an HSP: its range on the plus strand of the query
/// (`begin`, already folded for a minus context by the caller, and `len`, the query extent of
/// the HSP), its raw score and its e-value.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:233-243
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct BestHitKey {
    pub(crate) begin: i32,
    pub(crate) len: i32,
    pub(crate) score: i32,
    pub(crate) evalue: f64,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:42-49
/// `hsp` is the slot of the payload in the arena. `begin` and `end` are the range with the
/// overhang of the stage (or without it after the list was rebuilt, `s_ImportFromHitlist`);
/// `key_begin` is the range's start without it, which a rebuild reads again from the HSP.
#[derive(Clone, Copy)]
struct LinkedHspBh {
    hsp: usize,
    sid: u32,
    begin: i64,
    end: i64,
    len: i32,
    key_begin: i32,
    score: i32,
    evalue: f64,
}

/// The HSPs that the filter keeps for one query: one list per subject, in the order of the
/// subject's first HSP in the query's list, each list sorted by score.
pub(crate) type BestHitQuery<H> = Vec<(u32, Vec<H>)>;

/// NCBI's best hit writer (or pipe) over the HSPs of one search: `s_BlastHSPBestHitInit`
/// (`writer`, `pipe`), `s_BlastHSPBestHitRun` (`save`, once per HSP) and
/// `s_BlastHSPBestHitFinal` (`finalize`).
///
/// NCBI frees a dropped or removed HSP at once; the arena reuses its slot for the next HSP, so it
/// holds at most as many HSPs as were kept together.
pub(crate) struct BestHitLists<H> {
    overhang: f64,
    /// `param_s`, one minus the score edge.
    param_s: f64,
    /// The list of each query (`best_list`), by query index.
    lists: Vec<Vec<LinkedHspBh>>,
    /// `num_hsps` and `max_hsps` of each query.
    counts: Vec<(i64, i64)>,
    /// The first `max_hsps` of a query: twice `prelim_hitlist_size`.
    first_max: i64,
    hsps: Vec<Option<H>>,
    free_hsps: Vec<usize>,
}

impl<H> BestHitLists<H> {
    /// The writer of the preliminary stage: the overhang is 0.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329
    pub(crate) fn writer(params: BestHitParams, prelim_hitlist_size: usize) -> Self {
        Self::new(0.0, params.score_edge, prelim_hitlist_size)
    }

    /// The pipe after the traceback: both values.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:377-385
    pub(crate) fn pipe(params: BestHitParams, prelim_hitlist_size: usize) -> Self {
        Self::new(params.overhang, params.score_edge, prelim_hitlist_size)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:65-78
    fn new(overhang: f64, score_edge: f64, prelim_hitlist_size: usize) -> Self {
        let size = i64::try_from(prelim_hitlist_size).unwrap_or(i64::MAX);
        Self {
            overhang,
            param_s: 1.0 - score_edge,
            lists: Vec::new(),
            counts: Vec::new(),
            first_max: size.saturating_mul(2),
            hsps: Vec::new(),
            free_hsps: Vec::new(),
        }
    }

    fn store(&mut self, hsp: H) -> usize {
        match self.free_hsps.pop() {
            Some(slot) => {
                self.hsps[slot] = Some(hsp);
                slot
            }
            None => {
                self.hsps.push(Some(hsp));
                self.hsps.len() - 1
            }
        }
    }

    fn free(&mut self, slot: usize) {
        self.hsps[slot] = None;
        self.free_hsps.push(slot);
    }

    /// `s_BlastHSPBestHitRun` for one HSP of subject `sid`, in the list of `query`. The HSPs
    /// must come in the order NCBI writes them (the subjects in OID order, each subject's HSPs
    /// in the order of its list), because the result depends on it.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:228-332
    pub(crate) fn save(&mut self, query: usize, sid: u32, key: BestHitKey, hsp: H) {
        if self.lists.len() <= query {
            self.lists.resize_with(query + 1, Vec::new);
            self.counts.resize(query + 1, (0, self.first_max));
        }
        let begin = i64::from(key.begin);
        let len_a = key.len;
        let end = begin + i64::from(len_a);
        let score_a = key.score;
        let evalue_a = key.evalue;
        let param_s = self.param_s;
        let den_a = f64::from(score_a) / f64::from(len_a) / param_s;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:245-262
        // Is the new HSP bad: does an HSP of the list hold it with a better e-value and score?
        let list = &self.lists[query];
        let mut at = 0;
        while at < list.len() && list[at].end < end {
            at += 1;
        }
        while at < list.len() && list[at].begin <= begin {
            let p = &list[at];
            if p.end >= end && p.evalue <= evalue_a && f64::from(p.score) / f64::from(p.len) > den_a
            {
                return;
            }
            at += 1;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:264-294
        // The HSPs of the list that the new HSP makes bad are removed.
        let param_overhang = self.overhang;
        let overhang = i64::from(ncbi_int4_from_double(
            2.0 * f64::from(len_a) * param_overhang / (1.0 - 2.0 * param_overhang),
        ));
        let allowed_begin = begin - overhang;
        let allowed_end = end + overhang;
        let overhang = i64::from(ncbi_int4_from_double(f64::from(len_a) * param_overhang));
        let begin = begin - overhang;
        let end = end + overhang;
        let den_a = f64::from(score_a) / f64::from(len_a) * param_s;
        let list = &mut self.lists[query];
        let start = list.partition_point(|p| p.begin < allowed_begin);
        let stop = list.partition_point(|p| p.begin < allowed_end);
        let bad = |p: &LinkedHspBh| {
            let overhang = (p.end - p.begin - i64::from(p.len)) / 2;
            p.begin + overhang >= begin
                && p.end - overhang <= end
                && p.evalue >= evalue_a
                && f64::from(p.score) / f64::from(p.len) < den_a
        };
        let removed: Vec<usize> = list[start..stop.max(start)]
            .iter()
            .filter(|p| bad(p))
            .map(|p| p.hsp)
            .collect();
        if !removed.is_empty() {
            let mut index = 0;
            list.retain(|p| {
                let keep = index < start || index >= stop || !bad(p);
                index += 1;
                keep
            });
            self.counts[query].0 -= removed.len() as i64;
            for slot in removed {
                self.free(slot);
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:296-320
        // The new HSP goes before the first HSP whose start is not smaller.
        let slot = self.store(hsp);
        let list = &mut self.lists[query];
        let at = list.partition_point(|p| p.begin < begin);
        list.insert(
            at,
            LinkedHspBh {
                hsp: slot,
                sid,
                begin,
                end,
                len: len_a,
                key_begin: key.begin,
                score: score_a,
                evalue: evalue_a,
            },
        );
        let (num_hsps, max_hsps) = &mut self.counts[query];
        *num_hsps += 1;
        if *num_hsps > *max_hsps {
            self.rebuild(query);
        }
    }

    /// `s_ExportToHitlist` then `s_ImportFromHitlist` on the list of `query`: the HSPs grouped by
    /// subject (in the order of each subject's first HSP) and inserted again by their start
    /// without the overhang, a later HSP before an earlier one of the same start. Nothing is
    /// dropped (an HSP list of NCBI's applications holds any number of HSPs, `BlastHspNumMax`),
    /// and the next rebuild waits for twice as many HSPs.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:85-182
    fn rebuild(&mut self, query: usize) {
        let list = std::mem::take(&mut self.lists[query]);
        let mut grouped: Vec<Vec<LinkedHspBh>> = Vec::new();
        let mut position: HashMap<u32, usize> = HashMap::new();
        for entry in list {
            match position.get(&entry.sid) {
                Some(&at) => grouped[at].push(entry),
                None => {
                    position.insert(entry.sid, grouped.len());
                    grouped.push(vec![entry]);
                }
            }
        }
        // Each insertion goes before the HSPs of the same start: the order of the insertions
        // reversed, then sorted (stably) by start.
        let mut rebuilt: Vec<LinkedHspBh> = grouped.into_iter().flatten().collect();
        rebuilt.reverse();
        rebuilt.sort_by_key(|entry| entry.key_begin);
        for entry in &mut rebuilt {
            entry.begin = i64::from(entry.key_begin);
            entry.end = entry.begin + i64::from(entry.len);
        }
        let count = rebuilt.len() as i64;
        self.lists[query] = rebuilt;
        self.counts[query] = (count, count.saturating_mul(2));
    }

    /// `s_BlastHSPBestHitFinal` up to its hit list: for each query with a list, one list per
    /// subject in the order of the subject's first HSP in the query's list, each sorted by
    /// `cmp` (`Blast_HSPListSortByScore`: the program's `ScoreCompareHSPs`). The caller sorts
    /// the lists by e-value and saves them into a hit list of `prelim_hitlist_size`
    /// (`Blast_HitListSortByEvalue`, `Blast_HitListUpdate`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:190-221
    pub(crate) fn finalize(
        mut self,
        mut cmp: impl FnMut(&H, &H) -> Ordering,
    ) -> Vec<(usize, BestHitQuery<H>)> {
        let lists = std::mem::take(&mut self.lists);
        let mut out = Vec::new();
        for (query, list) in lists.into_iter().enumerate() {
            if list.is_empty() {
                continue;
            }
            let mut lists: BestHitQuery<H> = Vec::new();
            let mut position: HashMap<u32, usize> = HashMap::new();
            for entry in list {
                let hsp = self.hsps[entry.hsp]
                    .take()
                    .expect("an HSP is in one list once");
                match position.get(&entry.sid) {
                    Some(&at) => lists[at].1.push(hsp),
                    None => {
                        position.insert(entry.sid, lists.len());
                        lists.push((entry.sid, vec![hsp]));
                    }
                }
            }
            for (_, hsps) in &mut lists {
                hsps.sort_by(&mut cmp);
            }
            out.push((query, lists));
        }
        out
    }
}

#[cfg(test)]
#[path = "hspfilter_besthit_tests.rs"]
mod tests;
