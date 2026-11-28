//! Internal BLASTX HSP ownership, collector, stream and result reduction.
use super::{
    args::ResolvedOptions,
    linking::LinkedHspList,
    preliminary::{compare_score, PreliminaryHsp},
    split::QueryChunk,
};
use crate::common::GapEditOp;
use anyhow::{ensure, Result};
use std::cmp::Ordering;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3298
pub trait ResultsObserver {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:878-894
    fn targets(&mut self, _stage: &str, _q: usize, _lists: &[HspList]) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676
    fn reap(&mut self, _stage: &str, _oid: usize, _list: &LinkedHspList, _cutoff: f64) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:399-433
    fn merge(&mut self, _chunk: &QueryChunk, _enter: bool) {}
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:268-281
    fn read(&mut self, _list: Option<&HspList>) {}
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:591-615
    fn batch_read(&mut self, _lists: &[HspList]) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:931-940
    fn containment(&mut self, _h: &PreliminaryHsp, _contained: bool) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-568
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-150
    fn cull_compare(&mut self, _ao: usize, _a: &Hsp, _bo: usize, _b: &Hsp, _result: bool) {}
    fn cull_delete(&mut self, _oid: usize, _hsp: &Hsp, _merit: i64) {}
    fn cull_save(&mut self, _oid: usize, _hsp: &Hsp) {}
    fn cull_saved(&mut self, _result: bool) {}
    fn cull_fork(&mut self, _begin: i32, _end: i32) {}
    fn cull_final(&mut self, _num_queries: usize) {}
    fn preliminary_call(&mut self, _oid: usize) {}
    fn numeric(&mut self, _stage: &str, _oid: usize, _hsps: &[PreliminaryHsp]) {}
    fn hsps(&mut self, _stage: &str, _oid: usize, _hsps: &[Hsp]) {}
    fn hitlist(
        &mut self,
        _stage: &str,
        _oid: usize,
        _max: usize,
        _lists: &[HspList],
        _candidate: Option<&HspList>,
    ) {
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1917-1934
#[derive(Clone, Debug)]
pub struct Hsp {
    pub hsp: PreliminaryHsp,
    pub num: i32,
    pub evalue: f64,
    pub bit_score: f64,
    pub identity: usize,
    pub positive: usize,
    pub edit_script: Vec<GapEditOp>,
    pub composition_method: i32,
    /// The query letters that `-subject`'s Seq-loc mapper left out of the HSP
    /// (`HspList::map_subject_alignments`: an end `Ins` segment), which its `qcovhsp` and its
    /// Seq-align's `hsp_percent_coverage` still count
    /// (`serial::seqalign::MappedScript::query_left_out`); 0 for every other HSP.
    pub query_left_out: u32,
}

impl Hsp {
    /// The ends of the HSP's Seq-align, which every report prints: the 0-based starts plus the
    /// letters of the edit script (`Sub` and `Ins` on the query, `Sub` and `Del` on the
    /// subject), or the HSP's ends without one. They are the HSP's own ends except where NCBI's
    /// gapped traceback had no right extension (`found_end` FALSE: rpstblastn's window copy of a
    /// frame shorter than the HSP's gapped start, `-strand plus` or `minus` with
    /// composition-based statistics). There the HSP keeps NCBI's ends for every step on the HSP,
    /// the sort before the Seq-aligns (`compare_hsp`) among them.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4650-4653
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:181-182 (the
    /// Seq-align runs the edit script from the starts)
    pub(crate) fn seqalign_ends(&self) -> (i32, i32) {
        if self.edit_script.is_empty() {
            return (self.hsp.q_end, self.hsp.s_end);
        }
        let (mut q, mut s) = (self.hsp.q_start, self.hsp.s_start);
        for op in &self.edit_script {
            match *op {
                GapEditOp::Sub(n) => {
                    q += n as i32;
                    s += n as i32;
                }
                GapEditOp::Ins(n) => q += n as i32,
                GapEditOp::Del(n) => s += n as i32,
            }
        }
        (q, s)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1390-1435
pub(crate) fn compare_evalue(a: f64, b: f64) -> Ordering {
    if a < 1.0e-180 && b < 1.0e-180 {
        Ordering::Equal
    } else {
        a.partial_cmp(&b).expect("finite BLAST E-value")
    }
}
pub(crate) fn compare_hsp(a: &Hsp, b: &Hsp) -> Ordering {
    compare_evalue(a.evalue, b.evalue).then_with(|| compare_score(&a.hsp, &b.hsp))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3078-3106
#[derive(Clone, Debug)]
pub struct HspList {
    pub query_index: usize,
    pub oid: usize,
    pub hsps: Vec<Hsp>,
    pub best_evalue: f64,
}
impl HspList {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1737-1749
    pub(crate) fn refresh(&mut self) {
        self.best_evalue = crate::common::best_evalue(self.hsps.iter().map(|h| h.evalue));
    }
    pub(crate) fn sort_evalue(&mut self) {
        self.hsps.sort_by(compare_hsp);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1632-1638
    /// `-subject`'s Seq-loc mapper (`s_RemapToSubjectLoc`) on the list's HSPs, after the report's
    /// sort (`sort_evalue`; the mapper works on the sorted Seq-aligns): each HSP as
    /// `serial::seqalign::mapped_script` keeps it, an end segment with a gap row left out (a
    /// Smith-Waterman alignment of `-use_sw_tback` that starts with a gap): its starts move past
    /// the left-out letters and its ends back before them; the scores and the identities stay, and
    /// so does the coverage (`Hsp::query_left_out`). Every format reads the mapped HSP. Only
    /// `BatchResults::put_in_report_order` maps, once per batch: a second mapping would leave out
    /// a segment again (a script that starts `[Del, Ins, ..]`).
    fn map_subject_alignments(&mut self) {
        for hsp in &mut self.hsps {
            let mapped = crate::report::serial::seqalign::mapped_script(&hsp.edit_script);
            if !mapped.drops(&hsp.edit_script) {
                continue;
            }
            let raw = &mut hsp.hsp;
            raw.q_start += mapped.start[0] as i32;
            raw.s_start += mapped.start[1] as i32;
            raw.q_end -= mapped.end[0] as i32;
            raw.s_end -= mapped.end[1] as i32;
            hsp.query_left_out = mapped.query_left_out();
            hsp.edit_script = hsp.edit_script[mapped.kept].to_vec();
        }
    }

    /// The HSPs in the order `sort_evalue` puts them (the reports' order, `report::sort_lists`),
    /// without moving them: the same stable sort by `compare_hsp`.
    pub(crate) fn evalue_order(&self) -> Vec<&Hsp> {
        let mut hsps: Vec<&Hsp> = self.hsps.iter().collect();
        hsps.sort_by(|a, b| compare_hsp(a, b));
        hsps
    }
}
impl super::runtime::BatchResults {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
    /// Puts the batch's lists in the order of every report, once: each list's HSPs sorted by
    /// e-value (NCBI sorts each list's HSPs, on their own ends, before it makes their Seq-aligns;
    /// `Hsp::seqalign_ends`), then, where the subjects are `-subject`'s records
    /// (`SubjectLoc::Interval`), mapped as its Seq-loc mapper maps them
    /// (`HspList::map_subject_alignments`). A later call leaves the lists as they are
    /// (`in_report_order`): the mapper leaves out a first and a last segment with a gap row once
    /// each, and a sort of the mapped HSPs could order ties otherwise.
    pub(crate) fn put_in_report_order(&mut self, loc: super::empty_rows::SubjectLoc) {
        if std::mem::replace(&mut self.in_report_order, true) {
            return;
        }
        for lists in &mut self.queries {
            for list in lists {
                list.sort_evalue();
                if loc == super::empty_rows::SubjectLoc::Interval {
                    list.map_subject_alignments();
                }
            }
        }
    }
}

pub(crate) fn compare_list(a: &HspList, b: &HspList) -> Ordering {
    match (a.hsps.first(), b.hsps.first()) {
        (None, None) => Ordering::Equal,
        (None, Some(_)) => Ordering::Greater,
        (Some(_), None) => Ordering::Less,
        (Some(ah), Some(bh)) => compare_evalue(a.best_evalue, b.best_evalue)
            .then_with(|| bh.hsp.score.cmp(&ah.hsp.score))
            .then_with(|| b.oid.cmp(&a.oid)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3298
pub(crate) struct HitList {
    pub lists: Vec<HspList>,
    max: usize,
    heapified: bool,
}
impl HitList {
    pub fn new(max: usize) -> Self {
        Self {
            lists: Vec::new(),
            max,
            heapified: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3328-3337
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:827-857
    // (`core::blast_hspstream`)
    /// `BlastHSPCBSStreamClose`'s trim of this query's preliminary hit list (`ADAPTIVE_CBS`):
    /// sorted by E-value (the order is total, `compare_list`; the lists left without HSPs, last,
    /// removed) and cut where the E-values pass the limit, when it holds enough lists.
    fn cbs_close_trim(&mut self, hitlist_size: usize) {
        use crate::core::blast_hspstream::{cbs_trim_count, cbs_trim_sorts};
        let count = self.lists.len();
        if !cbs_trim_sorts(count, hitlist_size) {
            return;
        }
        self.lists.sort_unstable_by(compare_list);
        let sorted = self
            .lists
            .iter()
            .position(|list| list.hsps.is_empty())
            .unwrap_or(count);
        self.lists.truncate(sorted);
        self.heapified = false;
        let lists = &self.lists;
        if let Some(keep) = cbs_trim_count(count, sorted, hitlist_size, |i| lists[i].best_evalue) {
            self.lists.truncate(keep);
        }
    }
    fn down(&mut self, mut i: usize) {
        while i * 2 + 1 < self.lists.len() {
            let l = i * 2 + 1;
            let r = l + 1;
            let worse = if r == self.lists.len()
                || compare_list(&self.lists[l], &self.lists[r]) != Ordering::Less
            {
                l
            } else {
                r
            };
            if compare_list(&self.lists[i], &self.lists[worse]) != Ordering::Less {
                break;
            }
            self.lists.swap(i, worse);
            i = worse;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3298
    pub fn update_observed(
        &mut self,
        candidate: HspList,
        trace: &mut dyn ResultsObserver,
    ) -> Option<HspList> {
        let oid = candidate.oid;
        trace.hitlist("IN", oid, self.max, &self.lists, Some(&candidate));
        let discarded = self.update(candidate);
        trace.hitlist("OUT", oid, self.max, &self.lists, None);
        discarded
    }
    pub fn update(&mut self, mut candidate: HspList) -> Option<HspList> {
        candidate.refresh();
        if self.lists.len() < self.max {
            self.lists.push(candidate);
            return None;
        }
        if !self.heapified {
            for l in &mut self.lists {
                l.sort_evalue();
            }
            for i in (0..self.lists.len() / 2).rev() {
                self.down(i);
            }
            self.heapified = true;
        }
        candidate.sort_evalue();
        if compare_list(&self.lists[0], &candidate) == Ordering::Less {
            return Some(candidate);
        }
        let discarded = std::mem::replace(&mut self.lists[0], candidate);
        self.down(0);
        Some(discarded)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:44-70
/// The size of `blastn/hsp.rs` `get_prelim_hitlist_size` (with `ADAPTIVE_CBS`), computed
/// without NCBI's `int` overflow.
pub(crate) fn prelim_hitlist_size(o: &ResolvedOptions) -> usize {
    crate::algorithm::blastn::hsp::get_prelim_hitlist_size(
        o.hitlist_size as usize,
        o.composition > 0,
        o.gapped,
    )
}

/// `-subject_besthit` over one subject's HSPs of BLASTX (`core::blast_hits::subject_besthit`):
/// the contexts are the query frames, and the offsets those of the HSP in its frame.
pub(crate) fn subject_besthit(hsps: &mut Vec<Hsp>) {
    crate::core::blast_hits::subject_besthit(hsps, false, |h| best_hit_extent(&h.hsp));
}

/// `subject_besthit` of a list of preliminary HSPs sorted by score: the list of a sequence with
/// `-db_hard_mask` after each of its chunks joins it (`search_db_masked.rs`; series NP, M4b).
pub(crate) fn subject_besthit_preliminary(hsps: &mut Vec<PreliminaryHsp>) {
    crate::core::blast_hits::subject_besthit(hsps, false, best_hit_extent);
}

/// What `-subject_besthit` reads of a BLASTX HSP.
fn best_hit_extent(hsp: &PreliminaryHsp) -> crate::core::blast_hits::BestHitExtent {
    crate::core::blast_hits::BestHitExtent {
        context: hsp.context as i64,
        query_frame: i32::from(hsp.frame),
        offset: i64::from(hsp.q_start),
        end: i64::from(hsp.q_end),
        query_length: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:399-421
pub(crate) struct HspStream {
    pub queries: Vec<HitList>,
    sorted: Option<Vec<HspList>>,
    finalized: bool,
    culling: Option<super::culling::Culling>,
    /// The best hit writer of `-best_hit_overhang` / `-best_hit_score_edge`, in place of the
    /// collector (`with_best_hit`).
    besthit: Option<super::besthit::BestHitWriter>,
    max: usize,
}
impl HspStream {
    pub fn new(num_queries: usize, max: usize) -> Self {
        Self {
            queries: (0..num_queries).map(|_| HitList::new(max)).collect(),
            sorted: None,
            finalized: false,
            culling: None,
            besthit: None,
            max,
        }
    }

    /// The stream of `num_queries` queries whose writer is the best hit writer `writer`
    /// (`CSetupFactory::CreateHspWriter` with a best hit score edge), of `max` lists per query.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329
    pub(crate) fn with_best_hit(
        writer: super::besthit::BestHitWriter,
        num_queries: usize,
        max: usize,
    ) -> Self {
        let mut stream = Self::new(num_queries, max);
        stream.besthit = Some(writer);
        stream
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:105-161
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:330-342
    // NCBI's sum is an `Int4` that overflows for a limit of 2147483645 or more; LOSAT adds in
    // `i64` (NCBI parity DECISIONS 17), as `core::hspfilter_culling::CullingForest::prelim_max`.
    pub fn with_culling(lengths: Vec<i32>, limit: i32, max: usize) -> Self {
        let mut stream = Self::new(lengths.len() / 6, max);
        let limit = i64::from(limit);
        stream.culling = Some(super::culling::Culling::new(
            lengths,
            if limit > 1 { limit + 3 } else { limit },
        ));
        stream
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:604-647
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:105-161
    pub fn write_observed(
        &mut self,
        oid: usize,
        hsps: Vec<Hsp>,
        trace: &mut dyn ResultsObserver,
    ) -> Result<Vec<HspList>> {
        ensure!(
            !self.finalized && self.sorted.is_none(),
            "BLASTX HSP stream is finalized"
        );
        if let Some(culling) = &mut self.culling {
            culling.write(oid, hsps, trace);
            return Ok(Vec::new());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:237-352
        if let Some(writer) = &mut self.besthit {
            writer.write(oid, hsps);
            return Ok(Vec::new());
        }
        let mut lists: Vec<Vec<Hsp>> = (0..self.queries.len()).map(|_| Vec::new()).collect();
        for h in hsps {
            lists[h.hsp.context / 6].push(h);
        }
        let mut discarded = Vec::new();
        for (query_index, hsps) in lists.into_iter().enumerate() {
            if hsps.is_empty() {
                continue;
            }
            let list = HspList {
                query_index,
                oid,
                hsps,
                best_evalue: 0.0,
            };
            if let Some(l) = self.queries[query_index].update_observed(list, trace) {
                discarded.push(l);
            }
        }
        Ok(discarded)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:99-130
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:105-130
    pub fn finalize_observed(&mut self, trace: &mut dyn ResultsObserver) {
        if self.finalized {
            return;
        }
        if let Some(culling) = &mut self.culling {
            self.queries = culling.finalize(self.queries.len(), self.max);
            trace.cull_final(self.queries.len());
            for q in &self.queries {
                for l in &q.lists {
                    trace.hsps("CULL_FINAL", l.oid, &l.hsps);
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:194-231
        if let Some(writer) = self.besthit.take() {
            self.queries = writer.finalize(self.queries.len(), self.max);
        }
        self.finalized = true;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:456-533
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:419-421
    pub fn merge_observed(
        &mut self,
        mut chunk_stream: Self,
        chunk: &QueryChunk,
        trace: &mut dyn ResultsObserver,
    ) -> Result<Vec<HspList>> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:399-421
        trace.merge(chunk, true);
        ensure!(self.sorted.is_none(), "cannot merge a read BLASTX stream");
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:419-420
        chunk_stream.finalize_observed(trace);
        self.finalize_observed(trace);
        let mut discarded = Vec::new();
        for (local_query, mut incoming_hitlist) in chunk_stream.queries.into_iter().enumerate() {
            if incoming_hitlist.lists.is_empty() {
                continue;
            }
            let global_query = chunk.absolute_contexts[local_query * 6..local_query * 6 + 6]
                .iter()
                .copied()
                .find(|c| *c >= 0)
                .expect("query with HSP has context") as usize
                / 6;
            let mut points = [-1; 6];
            for local in local_query * 6..local_query * 6 + 6 {
                if chunk.absolute_contexts[local] >= 0 {
                    points[chunk.absolute_contexts[local] as usize % 6] = chunk.corrections[local];
                }
            }
            for list in &mut incoming_hitlist.lists {
                list.query_index = global_query;
                for h in &mut list.hsps {
                    let local = h.hsp.context;
                    ensure!(
                        chunk.absolute_contexts[local] >= 0,
                        "HSP in invalid split context"
                    );
                    h.hsp.context = chunk.absolute_contexts[local] as usize;
                    let offset = chunk.corrections[local];
                    h.hsp.q_start += offset;
                    h.hsp.q_end += offset;
                    h.hsp.q_gapped_start += offset;
                    h.hsp.frame = [1, 2, 3, -1, -2, -3][h.hsp.context % 6];
                }
            }
            let target = &mut self.queries[global_query];
            if target.lists.is_empty() {
                *target = incoming_hitlist;
                continue;
            }
            incoming_hitlist.lists.sort_by_key(|l| l.oid);
            target.lists.sort_by_key(|l| l.oid);
            let mut combined = HitList::new(incoming_hitlist.max);
            let mut old = std::mem::take(&mut target.lists).into_iter().peekable();
            let mut new = incoming_hitlist.lists.into_iter().peekable();
            while old.peek().is_some() || new.peek().is_some() {
                let list = match (old.peek(), new.peek()) {
                    (Some(a), Some(b)) if a.oid == b.oid => {
                        let mut a = old.next().unwrap();
                        let b = new.next().unwrap();
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-2884
                        if points.iter().any(|p| *p > 0) {
                            trace.hsps("OVERLAP_OLD", a.oid, &a.hsps);
                            trace.hsps("OVERLAP_NEW", b.oid, &b.hsps);
                            merge_hsps(&mut a.hsps, b.hsps, &points, chunk.overlap);
                            trace.hsps("OVERLAP_OUT", a.oid, &a.hsps);
                        } else {
                            a.hsps.extend(b.hsps);
                            a.hsps.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
                        }
                        a
                    }
                    (Some(a), Some(b)) if a.oid > b.oid => new.next().unwrap(),
                    (Some(_), _) => old.next().unwrap(),
                    _ => new.next().unwrap(),
                };
                if let Some(l) = combined.update_observed(list, trace) {
                    discarded.push(l);
                }
            }
            *target = combined;
        }
        for q in &mut self.queries {
            for l in &mut q.lists {
                l.hsps.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:521-537
        trace.merge(chunk, false);
        Ok(discarded)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:173-202
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:145-164
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:816-864
    /// The close with `ADAPTIVE_CBS` (`core::blast_hspstream::adaptive_cbs_close`): each query's
    /// hit list trimmed after the writer's final step, then `close_observed`.
    pub fn close_cbs_observed(&mut self, trace: &mut dyn ResultsObserver, hitlist_size: usize) {
        if self.sorted.is_some() {
            return;
        }
        self.finalize_observed(trace);
        for query in &mut self.queries {
            query.cbs_close_trim(hitlist_size);
        }
        self.close_observed(trace);
    }

    /// The lists of the closed stream, in their order, as `cut` leaves them (`Runtime::batch_end`:
    /// without the subjects that the fetches of the redo or the traceback exclude).
    pub fn retain_closed(&mut self, cut: impl FnOnce(&mut Vec<HspList>)) {
        if let Some(sorted) = self.sorted.as_mut() {
            cut(sorted);
        }
    }

    pub fn close_observed(&mut self, trace: &mut dyn ResultsObserver) {
        if self.sorted.is_some() {
            return;
        }
        self.finalize_observed(trace);
        let mut lists = Vec::new();
        for (query_index, q) in self.queries.iter_mut().enumerate() {
            for mut l in std::mem::take(&mut q.lists) {
                l.query_index = query_index;
                lists.push(l);
            }
        }
        lists.sort_by(|a, b| b.oid.cmp(&a.oid));
        self.sorted = Some(lists);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:320-329
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:145-164
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:268-281
    pub fn read_observed(&mut self, trace: &mut dyn ResultsObserver) -> Option<HspList> {
        self.close_observed(trace);
        let list = self.sorted.as_mut().unwrap().pop();
        trace.read(list.as_ref());
        list
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:599-613
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:145-164
    pub fn batch_read_observed(&mut self, trace: &mut dyn ResultsObserver) -> Vec<HspList> {
        self.close_observed(trace);
        let lists = self.sorted.as_mut().unwrap();
        let Some(oid) = lists.last().map(|l| l.oid) else {
            trace.batch_read(&[]);
            return Vec::new();
        };
        let mut batch = Vec::new();
        while lists.last().is_some_and(|l| l.oid == oid) {
            batch.push(lists.pop().unwrap());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:606-615
        trace.batch_read(&batch);
        batch
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2911-3003
fn merge_hsps(combined: &mut Vec<Hsp>, mut incoming: Vec<Hsp>, points: &[i32; 6], overlap: i32) {
    let mut nold = 0;
    for i in 0..combined.len() {
        let h = &combined[i].hsp;
        let p = points[h.context % 6];
        if p >= 0
            && if h.frame >= 0 {
                h.q_end > p
            } else {
                h.q_start < p + overlap
            }
        {
            combined.swap(nold, i);
            nold += 1;
        }
    }
    let mut nnew = 0;
    for i in 0..incoming.len() {
        let h = &incoming[i].hsp;
        let p = points[h.context % 6];
        if p >= 0
            && if h.frame < 0 {
                h.q_end > p
            } else {
                h.q_start < p + overlap
            }
        {
            incoming.swap(nnew, i);
            nnew += 1;
        }
    }
    let mut removed = vec![false; incoming.len()];
    for a in &mut combined[..nold] {
        for (i, b) in incoming[..nnew].iter().enumerate() {
            if removed[i] || a.hsp.context != b.hsp.context {
                continue;
            }
            let ah = &a.hsp;
            let bh = &b.hsp;
            let (end, start) = if ah.frame >= 0 {
                (ah.q_end - ah.s_end, bh.q_start - bh.s_start)
            } else {
                (bh.q_end - bh.s_end, ah.q_start - ah.s_start)
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1488-1531
            // Split-query searches are always gapped (`split_query_cxx.cpp:60-61`), so every HSP
            // has positive query length and protein score bounded by length * matrix max
            // (< Int4::MAX); `merge_two_chunk_hsps` is therefore always defined.
            if (end - start).abs() < 10
                && crate::core::blast_hits::merge_two_chunk_hsps(&mut a.hsp, &b.hsp, true).expect(
                    "gapped split-query protein HSPs have positive length and score bounded by \
                     length * matrix max (< Int4::MAX)",
                )
            {
                removed[i] = true;
            }
        }
    }
    combined.extend(
        incoming
            .into_iter()
            .enumerate()
            .filter_map(|(i, h)| (!removed[i]).then_some(h)),
    );
    combined.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:89-101
pub struct HeapEntry {
    pub list: HspList,
    pub best_score: i32,
}
fn worse(a: &HeapEntry, b: &HeapEntry) -> bool {
    a.list
        .best_evalue
        .partial_cmp(&b.list.best_evalue)
        .unwrap()
        .then_with(|| b.best_score.cmp(&a.best_score))
        .then_with(|| b.list.oid.cmp(&a.list.oid))
        == Ordering::Greater
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:402-419
pub struct CompoHeap {
    pub entries: Vec<HeapEntry>,
    threshold: usize,
    cutoff: f64,
    pub worst: f64,
    pub heapified: bool,
}
impl CompoHeap {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:163-163
    pub fn new(threshold: usize) -> Self {
        Self {
            entries: Vec::new(),
            threshold,
            cutoff: 0.002,
            worst: 0.0,
            heapified: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:163-245
    fn down(&mut self, mut i: usize) {
        loop {
            let l = i * 2 + 1;
            let r = l + 1;
            let mut large = i;
            if l < self.entries.len() && worse(&self.entries[l], &self.entries[i]) {
                large = l;
            }
            if r < self.entries.len() && worse(&self.entries[r], &self.entries[large]) {
                large = r;
            }
            if large == i {
                break;
            }
            self.entries.swap(i, large);
            i = large;
        }
    }
    fn up(&mut self, mut i: usize) {
        while i > 0 {
            let p = (i - 1) / 2;
            if !worse(&self.entries[i], &self.entries[p]) {
                break;
            }
            self.entries.swap(i, p);
            i = p;
        }
    }
    fn heapify(&mut self) {
        if !self.heapified {
            self.heapified = true;
            for i in (0..self.entries.len() / 2).rev() {
                self.down(i);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
    pub fn would_insert(&mut self, e: f64, score: i32, oid: usize) -> bool {
        if self.entries.len() < self.threshold || e <= self.cutoff || e < self.worst {
            return true;
        }
        self.heapify();
        let root = &self.entries[0];
        root.list
            .best_evalue
            .partial_cmp(&e)
            .unwrap()
            .then_with(|| score.cmp(&root.best_score))
            .then_with(|| oid.cmp(&root.list.oid))
            == Ordering::Greater
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:330-391
    pub fn insert(&mut self, entry: HeapEntry) -> Option<HeapEntry> {
        if !self.heapified && self.entries.len() >= self.threshold {
            self.heapify();
        }
        if !self.heapified {
            self.worst = self.worst.max(entry.list.best_evalue);
            self.entries.push(entry);
            return None;
        }
        let discarded = if self.entries.len() < self.threshold
            || (entry.list.best_evalue <= self.cutoff && self.worst <= self.cutoff)
        {
            self.entries.push(entry);
            self.up(self.entries.len() - 1);
            None
        } else if worse(&self.entries[0], &entry) {
            let d = std::mem::replace(&mut self.entries[0], entry);
            self.down(0);
            Some(d)
        } else {
            self.down(0);
            Some(entry)
        };
        self.worst = self.entries[0].list.best_evalue;
        discarded
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:439-461
    pub fn pop(&mut self) -> Option<HeapEntry> {
        self.heapify();
        if self.entries.is_empty() {
            return None;
        }
        let entry = self.entries.swap_remove(0);
        if !self.entries.is_empty() {
            self.down(0);
        }
        Some(entry)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1580
    /// The heap as `common::database::redo_fetches` reads it at a match's turn, `before` telling
    /// whether a record's subject came before the match in the redo.
    pub fn redo_heap(
        &self,
        before: impl Fn(usize) -> bool,
    ) -> crate::algorithm::common::database::RedoHeap {
        crate::algorithm::common::database::RedoHeap {
            ecutoff: self.cutoff,
            heap_threshold: self.threshold,
            inserted_before: self
                .entries
                .iter()
                .filter(|entry| entry.list.best_evalue <= self.cutoff && before(entry.list.oid))
                .count(),
            ends_early: true,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1580
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:404-409
    /// NCBI's early termination test of a match whose best E-value is `e`: every heap of `heaps` is
    /// filled to its cutoff (at least `threshold` entries, the worst at or below `cutoff`) and `e`
    /// is above five times the cutoff.
    ///
    /// A heap that `insert` filled from empty (a threshold of 1 or more) is filled to its cutoff
    /// exactly when at least `threshold` of its entries have a best E-value at or below `cutoff`.
    /// `worst` is the largest best E-value of the entries (the list's maximum, then the heap's
    /// root). As a list the heap holds at most `threshold` entries; as a heap it grows past
    /// `threshold` only by an entry at or below the cutoff added to entries all at or below it, and
    /// otherwise replaces its worst entry by a better one or keeps its entries. So each insert at or
    /// below the cutoff adds one entry at or below it (pushed, or in place of a worst entry above
    /// it), and no insert takes one away. `x_early_ahead` counts so; `results_boundary_tests.rs`
    /// checks the equivalence over random inserts.
    pub fn early(e: f64, heaps: &[Self]) -> bool {
        heaps.iter().all(|h| {
            h.entries.len() >= h.threshold
                && h.worst <= h.cutoff
                && crate::common::evalue_not_at_most(e, 5.0 * h.cutoff)
        })
    }

    /// No NCBI counterpart: scheduling only. The guess of the composition redo that offers matches
    /// to be redone ahead (`blastx/runtime/kappa_batch.rs`, `XPending`): `early(e, heaps)` on the
    /// heaps as they will be once each heap `heaps[q]` has taken `ahead[q]` more entries at or below
    /// its cutoff (the matches offered ahead of the loop that `x_counts_ahead` counts). That is,
    /// every heap holds at least `threshold` entries at or below `cutoff` counting `ahead[q]` more
    /// (see `early`), and `e` is above five times the cutoff. With `ahead` all zero it equals
    /// `early`. It never decides a result: the loop's test in stream order is `early`.
    #[cfg(any(test, feature = "parallel"))]
    pub(crate) fn x_early_ahead(e: f64, heaps: &[Self], ahead: &[usize]) -> bool {
        heaps.iter().enumerate().all(|(q, h)| {
            let needed = h.threshold.saturating_sub(ahead[q]);
            // A heap whose worst entry is at or below the cutoff has all its entries there; one
            // above it holds `threshold` entries at most, so the count reads that many.
            crate::common::evalue_not_at_most(e, 5.0 * h.cutoff)
                && h.entries.len() >= needed
                && (h.worst <= h.cutoff
                    || h.entries
                        .iter()
                        .filter(|entry| entry.list.best_evalue <= h.cutoff)
                        .take(needed)
                        .count()
                        >= needed)
        })
    }

    /// No NCBI counterpart: scheduling only. Whether a match offered ahead whose best E-value is
    /// `e` counts in `x_early_ahead`'s `ahead` for this heap (its query's): `e` at or below the
    /// cutoff, an entry the early termination test counts if the redo keeps it there. A guess
    /// either way (the redo moves the E-value), which changes the schedule only.
    #[cfg(any(test, feature = "parallel"))]
    pub(crate) fn x_counts_ahead(&self, e: f64) -> bool {
        e <= self.cutoff
    }

    /// The heap test of a Smith-Waterman window in e-value mode (`-use_sw_tback`, line RG item
    /// X4b-2), the one owner of the heap a window's query context reads: NCBI tests the heap that
    /// the window's query context `context` (6 x query + frame) indexes in the batch's array of
    /// `heaps`, which holds one heap per query. Below the batch's query count that is the heap of
    /// the batch's query number `context` (the window's own query only for context 0), whose
    /// answer this is. At or past it NCBI reads its other memory; LOSAT tests an empty heap there,
    /// which takes the window (approved exception 23 of PD-LOSAT-NCBI-DEFECTS, the Owner's
    /// decision Q-HEAP=B of 2026-11-21). NCBI's run of the batch followed by five queries without
    /// a hit per query, whose heaps are empty and inside the array, gives the same result where it
    /// searches the batch as the input as is: no heap of the batch fills to its cutoff (the early
    /// termination) and no lone query is split.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1465-1474
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3282-3296
    pub fn would_insert_at(
        heaps: &mut [Self],
        context: usize,
        e: f64,
        score: i32,
        oid: usize,
    ) -> bool {
        match heaps.get_mut(context) {
            Some(heap) => heap.would_insert(e, score, oid),
            // An empty heap takes every window: it holds fewer entries than its threshold, the
            // hit list size, which is at least 1.
            None => true,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
#[cfg(test)]
#[path = "results_boundary_tests.rs"]
mod boundary_tests;

#[cfg(test)]
#[path = "results_negative_length_tests.rs"]
mod negative_length_tests;

#[cfg(test)]
#[path = "results_report_order_tests.rs"]
mod report_order_tests;
