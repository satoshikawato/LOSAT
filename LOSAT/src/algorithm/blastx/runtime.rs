//! Native serial BLASTX runtime returning actual final HSPs and ancillary state.
use super::{
    args::ResolvedOptions,
    input::FastaRecord,
    kappa::KappaState,
    parameters::{subject_parameters, subject_update, ContextParameters},
    preliminary::PreliminaryHsp,
    query_setup::PreparedQueryBatch,
    reader_batch::UnitPrep,
    results::{
        prelim_hitlist_size, subject_besthit, CompoHeap, HeapEntry, HitList, Hsp, HspList,
        HspStream,
    },
    search::{search_core, PreliminarySink},
    split::QueryChunk,
    statistics::{
        evaluate_gapped, identity_and_positive, normalize_scores, reap_final, reap_preliminary,
    },
    subjects::Subjects,
};
use crate::core::blast_hits::{reap_by_query_coverage, trim_by_max_hsps};
use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::Result;
// The redo of a batch's matches with composition-based statistics (a child module, so that it uses
// the private types of this file).
mod kappa_batch;
use kappa_batch::KappaBatch;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1777
/// The query offset, query end and context length of an HSP (`Blast_HSPQueryCoverageTest`):
/// the offsets are in the residues of the HSP's query frame, whose length is the context's.
fn query_extent(
    hsp: &super::preliminary::PreliminaryHsp,
    contexts: &[super::query_setup::QueryContext],
) -> (i64, i64, i64) {
    (
        i64::from(hsp.q_start),
        i64::from(hsp.q_end),
        contexts[hsp.context].length as i64,
    )
}

pub struct BatchResults {
    pub query_ordinal: usize,
    pub queries: Vec<Vec<HspList>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-101
    pub prepared: PreparedQueryBatch,
    pub parameters: Vec<ContextParameters>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181,202-210
    pub search_skipped: bool,
    /// The ordinal id of the first of the batch's fetches with `check_oid_exclusion` (its
    /// traceback's or redo's, by subject) that read a header with a UTF8String where a
    /// VisibleString is expected: NCBI's warning comes there, before its Seq-aligns' reads
    /// (`ReportedSubjects::read_titles`; of the units of a reader batch, the least).
    pub fetched_utf8: Option<u32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:177-181 (quoted
    // at `report::archive::OptionsHandle::pass`)
    /// The `EffectiveSearchSpace` option that the batch's search leaves in the options handle
    /// where it is not context 0's search space (`archive::batch_searched_space`): RPSTBLASTN's
    /// last search unit's (`rps_search::volumes::RpsUnits::request_space`). `None`: the search
    /// space of context 0.
    pub request_space: Option<i64>,
    /// RPSTBLASTN: the HSPs of each query's profiles that the alias files' OID range leaves out
    /// past the formatter's cut, which left the query's merged lists
    /// (`rps_search::volumes::report_ordinal_ids`): outfmt 7's count and the prefetch's limit
    /// still count them (`BatchResults::unlisted_hsps`). Empty where none left.
    pub unlisted: Vec<usize>,
    /// RPSTBLASTN over several volumes: the number of each query's merged lists that the
    /// formatter's cut to the hit list size keeps, which counts following lists whose Seq-align
    /// ids match as one subject (`database::ReportedSubjects::for_batch`,
    /// `common::database::joined_cut`). Empty where the reports take
    /// `ResolvedOptions::shown_lists` (`BatchResults::cut_lists`).
    pub shown: Vec<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:102-104
    /// The gapped K of each query's footer where it is not the table's (`kbp_gap` of the query's
    /// first valid context): an RPSTBLASTN search with composition-based statistics writes NCBI's
    /// `kbp_gap[0]` (series PR, item PS12b, `rps_search::engine::batch_gapped_k`); `None` in
    /// BLASTX's results.
    pub gapped_k: Option<Vec<f64>>,
}

impl BatchResults {
    /// The HSPs of the query `query` of the batch whose profiles left its merged lists past the
    /// formatter's cut (`unlisted`), 0 for every other search.
    pub fn unlisted_hsps(&self, query: usize) -> usize {
        self.unlisted.get(query).copied().unwrap_or(0)
    }

    /// The lists of the query `query` that the reports NCBI's formatter cuts write
    /// (`ResolvedOptions::shown_lists` says which): its first `shown` lists, else its first
    /// `options.shown_lists()`.
    pub(crate) fn shown_lists(
        &self,
        query: usize,
        options: &super::args::ResolvedOptions,
    ) -> &[HspList] {
        self.cut_lists(query, options.shown_lists())
    }

    /// The lists of the query `query` that the formatter's cut keeps (`shown`, else the first
    /// `cut`, `ResolvedOptions::shown_lists`): the one reading of that cut, which
    /// `shown_lists` and `database::ReportedSubjects::for_batch` take.
    pub(crate) fn cut_lists(&self, query: usize, cut: usize) -> &[HspList] {
        let lists = &self.queries[query];
        &lists[..lists
            .len()
            .min(self.shown.get(query).copied().unwrap_or(cut))]
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1539-1554
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3647
pub trait RuntimeObserver: super::results::ResultsObserver + Send {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
    fn restored(&mut self, _batch: &PreparedQueryBatch, _params: &[ContextParameters]) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
    #[allow(clippy::too_many_arguments)]
    fn heap(
        &mut self,
        _query: usize,
        _stage: &str,
        _evalue: f64,
        _score: i32,
        _oid: i64,
        _decision: i64,
        _heap: &CompoHeap,
    ) {
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:757-764
    fn link_setup(
        &mut self,
        _subject_length: i32,
        _db_length: i64,
        _state: crate::core::link_hsps::LinkCutoffs,
        _link: &super::linking::LinkParameters,
    ) {
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
    fn link_run(
        &mut self,
        _subject_length: i32,
        _gapped: bool,
        _link: &super::linking::LinkParameters,
    ) {
    }
    fn early(&mut self, _evalue: f64, _queries: usize, _decision: bool) {}
    // No NCBI counterpart: tells the search whether any observer callback records something; it
    // does not change any value NCBI computes.
    /// True when no callback records anything, so that work may be done out of order (and some of
    /// it twice).
    fn x_passive(&self) -> bool {
        false
    }
}
struct Noop;
impl super::results::ResultsObserver for Noop {}
impl RuntimeObserver for Noop {
    fn x_passive(&self) -> bool {
        true
    }
}

// No NCBI counterpart: the rule of one parallel BLASTX stage and the record of the side the stage
// took (`LOSAT_X_REACH`); it does not change any value NCBI computes.
/// Whether a parallel stage runs: `rule` holds. Records the side taken. The stages that do not
/// run (one thread, no pool, a small batch) run the serial reference code.
#[cfg(feature = "parallel")]
pub(crate) fn x_bx_par_rule(rule: impl FnOnce() -> bool) -> bool {
    let own = rule();
    crate::utils::xreach::x_reach!("LOSAT_X_BXPAR", own);
    own
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334,3499-3502
// One XKappaWorker is the Rust counterpart of the per-thread state NCBI keeps in NRrecord_tld,
// redo_align_params_tld and matrix_tld. It is used by one thread at a time.
// One pool slot's Kappa state.  SAFETY of `Send`: the only non-Send field is the `context` cell of
// the gapping parameters, which `redo_context_observed` fills for one synchronous call and restores
// before returning; the state is only reached through its mutex.
#[cfg(feature = "parallel")]
struct XKappaWorker(KappaState);
#[cfg(feature = "parallel")]
unsafe impl Send for XKappaWorker {}

// No NCBI counterpart (NCBI keeps one record per thread for the whole search): how many redo
// workers a batch end makes. EXPERIMENT (LOSAT_X_BXKAPPAFIT, in the strict set: on unless `=0`):
// no more than its matches, since one speculative batch then holds every match and each match
// takes the worker of its position; `=0`: one per pool thread. A worker's state does not change a
// match's result (`x_speculate_match`: a result that read a predecessor's matrix is redone in
// stream order).
#[cfg(feature = "parallel")]
fn x_kappa_workers(threads: usize, matches: usize) -> usize {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if !*ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXKAPPAFIT")) {
        return threads;
    }
    crate::utils::xreach::x_reach!("LOSAT_X_BXKAPPAFIT", matches < threads);
    threads.min(matches)
}

// No NCBI counterpart: a placement of the same value. Boxed because a result slot of a speculative
// batch holds it: inline, each slot of `batch_end`'s `speculated` took 3.1 KiB, one batch's 512
// slots at 16 threads 1.6 MiB, and a slot written at each 3.1 KiB touches every page of the block.
// That block stayed resident in the glibc arena of each thread that ends a batch.
/// The matrix a match leaves for the match after it.
type XMatrixLeft = Box<crate::core::composition_adjustment::adjust_scores::AdjustedProteinMatrix>;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1234
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3626
// x_speculate_match redoes one match with the same steps. Blast_AdjustScores runs for a protein
// subject only when hsp_index == 0 (line 1233), so a match can read the matrix its predecessor left
// behind. X_REDO_PROBE records whether that happened; if it did the result is dropped and the match
// is redone in stream order.
/// Redo one match without knowing the matrix its predecessor leaves behind.  `None` when the match
/// would have used that matrix (it is then redone in order); otherwise the finished list, its best
/// score, and the matrix this match leaves (`None` if it leaves the old one).
#[cfg(feature = "parallel")]
fn x_speculate_match(
    worker: &mut KappaState,
    list: &HspList,
    batch: &PreparedQueryBatch,
    full: &[ContextParameters],
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
    db_length: i64,
) -> Option<(HspList, i32, Option<XMatrixLeft>)> {
    use crate::core::composition_adjustment::redo_alignment::X_REDO_PROBE;
    worker.x_set_matrix(None);
    X_REDO_PROBE.with(|probe| probe.set(0));
    let encoded = subjects.read_searched(list.oid, options).ok()?;
    let input: Vec<_> = list.hsps.iter().map(|h| h.hsp.clone()).collect();
    let redone = worker.redo_list_observed(&input, batch, &encoded).ok()?;
    let probe = X_REDO_PROBE.with(|probe| probe.get());
    if probe & 2 != 0 {
        return None;
    }
    let (output, best_score) = finalize_redone(
        &redone,
        batch,
        full,
        options,
        &encoded,
        db_length,
        list.query_index,
        list.oid,
        &mut Noop,
    )
    .ok()?;
    let matrix = if probe & 1 != 0 {
        worker.x_take_matrix().map(Box::new)
    } else {
        None
    };
    Some((output, best_score, matrix))
}
pub(crate) fn search_internal<'a>(
    records: &[FastaRecord],
    subjects: impl Into<Subjects<'a>>,
    options: &ResolvedOptions,
) -> Result<Vec<BatchResults>> {
    search_internal_observed(records, subjects, options, &mut Noop)
}
// No NCBI counterpart: when the queued batches (2+ threads, several query batches) are searched (as
// a flow instead of waves); the reports, their order and the values in them are those of the waves.
/// EXPERIMENT (LOSAT_X_BXFLOW): a switch of the strict set, on unless `LOSAT_X_BXFLOW=0`
///  (acts with several query batches at 2+ threads).
#[cfg(feature = "parallel")]
pub(crate) fn x_bx_flow() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXFLOW"))
}

thread_local! {
    // No NCBI counterpart: flag for the batched search (a thread that runs one whole batch must not
    // start parallel work of its own); it does not change any value NCBI computes.
    // Set while this thread runs a search that must not start parallel work of its own.
    static X_INNER_SERIAL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

// No NCBI counterpart: reads the flag above; it does not change any value NCBI computes.
/// True inside `x_search_internal_serial`.
pub(crate) fn x_inner_serial() -> bool {
    X_INNER_SERIAL.with(std::cell::Cell::get)
}

// No NCBI counterpart: scope guard that sets and restores the flag above; it does not change any
// value NCBI computes.
#[cfg(feature = "parallel")]
struct XInnerSerial(bool);

#[cfg(feature = "parallel")]
impl XInnerSerial {
    fn enter() -> Self {
        Self(X_INNER_SERIAL.with(|flag| flag.replace(true)))
    }
}

#[cfg(feature = "parallel")]
impl Drop for XInnerSerial {
    fn drop(&mut self) {
        X_INNER_SERIAL.with(|flag| flag.set(self.0));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:256-265
// NCBI runs the single-threaded search (CPrelimSearchRunner) when IsMultiThreaded() is false. This
// is search_internal run on one thread, whatever thread that is, so that several batches can run
// side by side.
/// `search_internal` as the one-thread search, whatever thread it is called on. The caller runs
/// several of these side by side; each takes the serial path of every stage (the Web's, `web.rs`).
#[cfg(all(feature = "parallel", any(target_arch = "wasm32", test)))]
pub(crate) fn x_search_internal_serial<'a>(
    records: &[FastaRecord],
    subjects: impl Into<Subjects<'a>>,
    options: &ResolvedOptions,
) -> Result<Vec<BatchResults>> {
    x_search_unit_serial(records, subjects, options, None)
}

/// `x_search_internal_serial` for one unit of a reader batch of several, prepared in phase 1
/// (`reader_batch.rs`); `None`: `records` are the reader batch (or the whole batch is one unit).
/// Its callers are the several-batches-at-a-time paths, which need a pool (`parallel`).
#[cfg(feature = "parallel")]
pub(crate) fn x_search_unit_serial<'a>(
    records: &[FastaRecord],
    subjects: impl Into<Subjects<'a>>,
    options: &ResolvedOptions,
    prep: Option<UnitPrep>,
) -> Result<Vec<BatchResults>> {
    let _serial = XInnerSerial::enter();
    x_search_unit_in_pool(
        records,
        subjects,
        options,
        &crate::utils::threading::SearchPool::x_serial(),
        prep,
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// The search of one query batch (lcl_blast.Run()) on a pool that the caller keeps. The search
// itself is search_internal.
/// `search_internal` on a pool that the caller keeps for all its query batches.
pub(crate) fn x_search_internal_in_pool<'a>(
    records: &[FastaRecord],
    subjects: impl Into<Subjects<'a>>,
    options: &ResolvedOptions,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<Vec<BatchResults>> {
    x_search_unit_in_pool(records, subjects, options, pool, None)
}

/// `x_search_internal_in_pool` for one unit of a reader batch of several, prepared in phase 1
/// (`reader_batch.rs`); `None`: `records` are the reader batch (or the whole batch is one unit).
pub(crate) fn x_search_unit_in_pool<'a>(
    records: &[FastaRecord],
    subjects: impl Into<Subjects<'a>>,
    options: &ResolvedOptions,
    pool: &crate::utils::threading::SearchPool<'_>,
    prep: Option<UnitPrep>,
) -> Result<Vec<BatchResults>> {
    let subjects = subjects.into();
    let mut observer = Noop;
    let mut pipeline = Runtime {
        observer: &mut observer,
        ungapped_link_state: None,
        options,
        subjects,
        x_cursor: subjects.x_cursor(),
        x_residues: subjects.x_residue_pass(),
        // The Gumbel block's length of the e-values (`SubjectNumbers::gumbel_total`).
        db_length: subjects.numbers()?.gumbel_total,
        ordinal: 0,
        stream: None,
        chunk: None,
        output: Vec::new(),
    };
    search_core(records, subjects, options, &mut pipeline, pool, prep)?;
    Ok(pipeline.output)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3647
pub(crate) fn search_internal_observed<'a>(
    records: &[FastaRecord],
    subjects: impl Into<Subjects<'a>>,
    options: &ResolvedOptions,
    observer: &mut dyn RuntimeObserver,
) -> Result<Vec<BatchResults>> {
    let subjects = subjects.into();
    let mut pipeline = Runtime {
        observer,
        ungapped_link_state: None,
        options,
        subjects,
        x_cursor: subjects.x_cursor(),
        x_residues: subjects.x_residue_pass(),
        // The Gumbel block's length of the e-values (`SubjectNumbers::gumbel_total`).
        db_length: subjects.numbers()?.gumbel_total,
        ordinal: 0,
        stream: None,
        chunk: None,
        output: Vec::new(),
    };
    crate::utils::threading::with_search_pool(options.num_threads as usize, "BLASTX", |pool| {
        search_core(records, subjects, options, &mut pipeline, pool, None)
    })?;
    Ok(pipeline.output)
}
struct Runtime<'a> {
    ungapped_link_state: Option<crate::core::link_hsps::LinkCutoffs>,
    observer: &'a mut dyn RuntimeObserver,
    options: &'a ResolvedOptions,
    subjects: Subjects<'a>,
    /// The lengths that the sink reads of a database read in
    /// blocks.
    x_cursor: Option<super::x_db_stream::XCursor<'a>>,
    /// The residues that the `-ungapped` sink reads of a database
    /// read in blocks.
    x_residues: Option<super::x_db_stream::XPass<'a>>,
    db_length: i64,
    ordinal: usize,
    stream: Option<HspStream>,
    chunk: Option<HspStream>,
    output: Vec<BatchResults>,
}

impl Runtime<'_> {
    /// The length of subject `oid`, read in the sink's order (for a `-db`: from the
    /// index, `x_db_stream::XCursor`).
    fn subject_length(&mut self, oid: usize) -> Result<usize> {
        match &mut self.x_cursor {
            Some(cursor) => cursor.length(oid),
            None => Ok(self.subjects.length(oid)),
        }
    }

    /// `PreliminarySink::subject` for the statistics of each subject: the subject's own
    /// parameters (`parameters::subject_update`, at the subject's length) in place of the
    /// chunk's `params`.
    #[cold]
    #[inline(never)]
    fn subject_own(
        &mut self,
        oid: usize,
        input: &[PreliminaryHsp],
        batch: &PreparedQueryBatch,
        params: &[ContextParameters],
        option_spaces: &[i64],
    ) -> Result<()> {
        let length = self.subject_length(oid)?;
        let own = subject_update(
            batch,
            params,
            self.options,
            self.subjects.numbers()?,
            length,
            option_spaces,
        )?;
        self.subject_with(oid, input, batch, &own)
    }

    /// The body of `PreliminarySink::subject` with the parameters `params` of the subject.
    fn subject_with(
        &mut self,
        oid: usize,
        input: &[PreliminaryHsp],
        batch: &PreparedQueryBatch,
        params: &[ContextParameters],
    ) -> Result<()> {
        let mut hsps: Vec<_> = input
            .iter()
            .map(|h| Hsp {
                hsp: h.clone(),
                num: 0,
                evalue: 0.0,
                bit_score: 0.0,
                identity: 0,
                positive: 0,
                edit_script: Vec::new(),
                composition_method: 0,
            })
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2538-2588
        if self.options.subject_besthit && !hsps.is_empty() {
            self.observer.hsps("SUBJECT_BEST_IN", oid, &hsps);
            subject_besthit(&mut hsps);
            self.observer.hsps("SUBJECT_BEST_OUT", oid, &hsps);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899
        let input: Vec<_> = hsps.iter().map(|h| h.hsp.clone()).collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1490-1542
        if !self.options.gapped {
            // With the database read in blocks, the length from
            // the index (`XCursor`), and the residues, which only `identity_and_positive` reads,
            // only of a subject with HSPs, from the sink's own pass (`XPass`).
            let x_length = if self.x_cursor.is_some() {
                Some(self.subject_length(oid)?)
            } else {
                None
            };
            let encoded = match &mut self.x_residues {
                Some(pass) if !input.is_empty() => pass.view(oid)?.encoded(oid),
                Some(_) => std::borrow::Cow::Borrowed(&[][..]),
                None => self.subjects.encoded(oid),
            };
            let subject_length = x_length.unwrap_or(encoded.len());
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
            if !input.is_empty() {
                if let Some(link) = super::linking::link_parameters(self.options) {
                    self.observer.link_run(subject_length as i32, false, &link);
                }
            }
            let mut linked = super::statistics::evaluate_ungapped(
                &input,
                batch,
                params,
                self.options,
                subject_length as i32,
                self.ungapped_link_state,
            )?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676
            let reap_cutoff = self.options.evalue
                * if self.options.composition > 1 && !self.options.old_fsc {
                    5.0
                } else {
                    1.0
                };
            self.observer.reap("PRE_REAP_IN", oid, &linked, reap_cutoff);
            reap_preliminary(&mut linked, self.options);
            self.observer
                .reap("PRE_REAP_OUT", oid, &linked, reap_cutoff);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2708-2736
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:776-790
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1448
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:74
            // BLASTX has protein subjects: ambiguity rescoring is inactive. Every HSP's query and
            // subject lengths are equal (`q_end - q_start`, NCBI's `Int4`), so its alignment length
            // is that difference. BLASTX shows no `-perc_identity` (`CGenericSearchArgs`' default
            // `show_perc_identity` false) and no minimal hit length (calloc's 0), so
            // `Blast_HSPTest` (`core::blast_hits::hsp_test`) runs with both at 0: it deletes the
            // HSPs of negative length that the one-hit word finder saves with `-window_size 0` and
            // a low `-threshold` (aa_ungapped.c:1054,1083), whose identity loop does not run, here,
            // after they took part in the first link and the reap above; it keeps those of length
            // 0, which are linked again and reported.
            let mut payloads: Vec<_> = linked
                .hsps
                .into_iter()
                .filter_map(|h| {
                    let align_length = h.hsp.q_end - h.hsp.q_start;
                    let (identity, positive) = match u32::try_from(align_length) {
                        Ok(length) => identity_and_positive(
                            &super::kappa::RedoneHsp {
                                hsp: h.hsp.clone(),
                                edit_script: vec![crate::common::GapEditOp::Sub(length)],
                                composition_method: 0,
                                identity_matrix: None,
                            },
                            batch,
                            &encoded,
                            self.options.scoring,
                        ),
                        Err(_) => (0, 0),
                    };
                    if crate::core::blast_hits::hsp_test(
                        identity as i64,
                        i64::from(align_length),
                        0.0,
                        0,
                    ) {
                        return None;
                    }
                    Some(Hsp {
                        hsp: h.hsp,
                        num: h.num,
                        evalue: h.evalue,
                        bit_score: 0.0,
                        identity,
                        positive,
                        edit_script: Vec::new(),
                        composition_method: 0,
                    })
                })
                .collect();
            payloads.sort_by(|a, b| super::preliminary::compare_score(&a.hsp, &b.hsp));
            let input: Vec<_> = payloads.iter().map(|h| h.hsp.clone()).collect();
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
            if !input.is_empty() {
                if let Some(link) = super::linking::link_parameters(self.options) {
                    self.observer.link_run(subject_length as i32, false, &link);
                }
            }
            let mut linked = super::statistics::evaluate_ungapped(
                &input,
                batch,
                params,
                self.options,
                subject_length as i32,
                self.ungapped_link_state,
            )?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676
            let reap_cutoff = self.options.evalue
                * if self.options.composition > 1 {
                    5.0
                } else {
                    1.0
                };
            self.observer.reap("PRE_REAP_IN", oid, &linked, reap_cutoff);
            reap_preliminary(&mut linked, self.options);
            self.observer
                .reap("PRE_REAP_OUT", oid, &linked, reap_cutoff);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1521-1538
            // An ungapped search removes the HSPs below the query coverage before the hit list
            // receives the subject (a subject left without HSPs never reaches it); the
            // traceback's reap then removes nothing, because NCBI never splits the query of an
            // ungapped search (split_query_cxx.cpp:60-61, "No split for ungapped mode"), so
            // both reaps read the same context lengths.
            if self.options.query_cov_hsp_perc > 0.0 {
                reap_by_query_coverage(&mut linked.hsps, self.options.query_cov_hsp_perc, |h| {
                    query_extent(&h.hsp, &batch.contexts)
                });
            }
            let hsps = linked
                .hsps
                .into_iter()
                .map(|h| {
                    let p = &params[h.hsp.context];
                    let mut payload = payloads[h.source_index].clone();
                    payload.bit_score = (h.hsp.score as f64 * p.ungapped.lambda
                        - p.ungapped.k.glibc_ln())
                        / std::f64::consts::LN_2;
                    payload.hsp = h.hsp;
                    payload.num = h.num;
                    payload.evalue = h.evalue;
                    payload
                })
                .collect();
            self.chunk
                .as_mut()
                .expect("chunk stream")
                .write_observed(oid, hsps, self.observer)?;
            return Ok(());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
        let subject_length = self.subject_length(oid)? as i32;
        if !input.is_empty() {
            if let Some(link) = super::linking::link_parameters(self.options) {
                self.observer.link_run(subject_length, true, &link);
            }
        }
        let mut linked = evaluate_gapped(
            &input,
            batch,
            params,
            self.options,
            subject_length,
            self.db_length,
            1.0,
        )?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676
        let reap_cutoff = self.options.evalue
            * if self.options.composition > 1 {
                5.0
            } else {
                1.0
            };
        self.observer.reap("PRE_REAP_IN", oid, &linked, reap_cutoff);
        reap_preliminary(&mut linked, self.options);
        self.observer
            .reap("PRE_REAP_OUT", oid, &linked, reap_cutoff);
        let hsps = linked
            .hsps
            .into_iter()
            .map(|h| {
                let mut payload = hsps[h.source_index].clone();
                payload.num = h.num;
                payload.evalue = h.evalue;
                payload
            })
            .collect();
        self.chunk
            .as_mut()
            .expect("chunk stream")
            .write_observed(oid, hsps, self.observer)?;
        Ok(())
    }
}
impl PreliminarySink for Runtime<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181,202-210
    fn unsearched_batch(
        &mut self,
        ordinal: usize,
        batch: &PreparedQueryBatch,
        params: &[ContextParameters],
    ) -> Result<()> {
        self.output.push(BatchResults {
            query_ordinal: ordinal,
            queries: vec![Vec::new(); batch.original_lengths.len()],
            prepared: batch.x_for_results(),
            parameters: params.to_vec(),
            search_skipped: true,
            fetched_utf8: None,
            request_space: None,
            unlisted: Vec::new(),
            shown: Vec::new(),
            gapped_k: None,
        });
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-235
    fn batch_start(
        &mut self,
        ordinal: usize,
        batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        self.ordinal = ordinal;
        self.stream = Some(if self.options.culling_limit > 0 {
            HspStream::with_culling(
                batch.contexts.iter().map(|c| c.length as i32).collect(),
                self.options.culling_limit,
                prelim_hitlist_size(self.options),
            )
        } else {
            HspStream::new(
                batch.original_lengths.len(),
                prelim_hitlist_size(self.options),
            )
        });
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-244
    fn chunk_start(&mut self, chunk: &QueryChunk) -> Result<()> {
        self.chunk = Some(if self.options.culling_limit > 0 {
            HspStream::with_culling(
                chunk
                    .prepared
                    .contexts
                    .iter()
                    .map(|c| c.length as i32)
                    .collect(),
                self.options.culling_limit,
                prelim_hitlist_size(self.options),
            )
        } else {
            HspStream::new(
                chunk.prepared.original_lengths.len(),
                prelim_hitlist_size(self.options),
            )
        });
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:583-593
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1450-1458
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-552
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3916-3921
    fn containment(&mut self, h: &PreliminaryHsp, contained: bool) {
        self.observer.containment(h, contained);
    }
    // No NCBI counterpart: tells the search whether the containment callback is read; it does not
    // change any value NCBI computes.
    fn x_wants_containment(&self) -> bool {
        !self.observer.x_passive()
    }
    fn preliminary_purge(
        &mut self,
        oid: usize,
        before: &[PreliminaryHsp],
        after: &[PreliminaryHsp],
    ) -> Result<()> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-568
        self.observer.preliminary_call(oid);
        // The score-only list is newly allocated with oid=0; the engine sets
        // its real subject OID only after endpoint purge/score sort.
        self.observer.numeric("PRE_ENDPOINT_IN", 0, before);
        self.observer.numeric("PRE_ENDPOINT_OUT", 0, after);
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455
    // search::compute_subject has performed this exact preparation before
    // WordFinder with worker-local state. Adopt it without recomputation,
    // retaining serial observer and linking state in input OID order.
    fn subject_start(
        &mut self,
        oid: usize,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
        prepared: Option<crate::core::link_hsps::LinkCutoffs>,
    ) -> Result<()> {
        self.ungapped_link_state = prepared;
        if let Some(prepared) = prepared {
            let length = self.subject_length(oid)? as i32;
            let link = super::linking::link_parameters(self.options)
                .expect("ungapped sum-stat parameters");
            self.observer
                .link_setup(length, self.db_length, prepared, &link);
        }
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1443 (cited above
    // `search::search_preliminary`): where NCBI computes the statistics of each subject (the
    // batch's `subject_option_spaces`), the e-values of a subject's HSPs read the subject's own
    // parameters, the ones `search::compute_own_subject` searched it with (`parameters::
    // subject_update` of the same chunk, options' spaces and length); a subject without HSPs
    // reads none. One test per subject on the default path.
    fn subject(
        &mut self,
        oid: usize,
        input: &[PreliminaryHsp],
        batch: &PreparedQueryBatch,
        params: &[ContextParameters],
    ) -> Result<()> {
        if let (Some(option_spaces), false) =
            (batch.subject_option_spaces.as_deref(), input.is_empty())
        {
            return self.subject_own(oid, input, batch, params, option_spaces);
        }
        self.subject_with(oid, input, batch, params)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:268-271
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:307-319
    fn chunk_end(&mut self, chunk: &QueryChunk, is_split: bool) -> Result<()> {
        let incoming = self.chunk.take().expect("chunk stream");
        if is_split {
            self.stream.as_mut().expect("batch stream").merge_observed(
                incoming,
                chunk,
                self.observer,
            )?;
        } else {
            self.stream = Some(incoming);
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
    fn batch_end(
        &mut self,
        batch: &PreparedQueryBatch,
        params: &[ContextParameters],
    ) -> Result<()> {
        // Diagnostic (LOSAT_TIMING): the clocks of the intervals of this function
        // (x_timing.rs). No NCBI counterpart.
        use super::x_timing::End as XEnd;
        let x_t = super::x_timing::lap();
        let mut x_t_finish = None;
        let mut stream = self.stream.take().expect("batch stream");

        let mut full = params.to_vec();
        subject_parameters(
            batch,
            &mut full,
            self.options,
            self.subjects.numbers()?,
            None,
        )?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1839-1846
        let hitlist_size = self.options.hitlist_size as usize;
        if crate::core::blast_hspstream::adaptive_cbs_close(
            self.options.composition > 0,
            hitlist_size,
        ) {
            stream.close_cbs_observed(self.observer, hitlist_size);
        } else {
            stream.close_observed(self.observer);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3560-3575 and
        // blast_traceback.c:1567-1577 (cited at `blastp::blast_engine::database::exclude_fetched`)
        // The redo's and the traceback's fetches (`check_oid_exclusion`) exclude a subject that the
        // negative seqid list keeps without a line (`FetchExcluded`): its lists leave the closed
        // stream before the redo or the traceback reads it. An ungapped search fetches no protein
        // subject; its such hits leave at the Seq-aligns (`ReportedSubjects::drop_unaligned`).
        // A header that the fetch cannot read ends NCBI's search: LOSAT's refusal
        // (`FetchExcluded::refusal`).
        // The subjects whose header holds a UTF8String where a VisibleString is expected, with
        // the query of their list, its best E-value and whether it has HSPs, and the lists left in
        // the stream's order (their query and subject) where there is one: the batch's first such
        // fetch is NCBI's warning (`fetched_utf8`; the redo's, `kappa_batch::redo_fetched`).
        let mut fetched_reads = Vec::new();
        let mut fetched_order = Vec::new();
        let mut fetched_utf8 = None;
        if self.options.gapped {
            if let Some(excluded) = self
                .subjects
                .seqdb()
                .and_then(crate::algorithm::common::database::FetchExcluded::of)
            {
                stream
                    .retain_closed(|lists| {
                        let found = excluded.retain_noting(
                            lists,
                            |list| u32::try_from(list.oid).ok(),
                            |list| (list.query_index, list.best_evalue, !list.hsps.is_empty()),
                        )?;
                        if !found.is_empty() {
                            fetched_order = lists
                                .iter()
                                .map(|list| (list.query_index, list.oid))
                                .collect();
                        }
                        fetched_reads = found;
                        Ok(())
                    })
                    .map_err(|failed| {
                        crate::algorithm::common::database::FetchExcluded::refusal(failed, "BLASTX")
                    })?;
                // The traceback fetches every subject of the stream.
                if self.options.composition == 0 {
                    fetched_utf8 = crate::algorithm::common::database::FetchedUtf8::traceback_first(
                        &fetched_reads
                            .iter()
                            .map(|(read, _)| *read)
                            .collect::<Vec<_>>(),
                    );
                }
            }
        }
        self.observer.restored(batch, &full);

        let mut queries: Vec<_> = (0..batch.original_lengths.len())
            .map(|_| HitList::new(self.options.hitlist_size as usize))
            .collect();

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3281-3296
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3120-3132
        let mut inactive_heaps: Vec<_> = if self.options.composition != 0 {
            (0..queries.len())
                .map(|_| CompoHeap::new(self.options.hitlist_size as usize))
                .collect()
        } else {
            Vec::new()
        };
        for (q, heap) in inactive_heaps.iter().enumerate() {
            self.observer.heap(q, "INIT", 0.0, 0, -1, 0, heap);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3281-3296
        let mut heaps: Vec<_> = if self.options.composition != 0 {
            (0..queries.len())
                .map(|_| CompoHeap::new(self.options.hitlist_size as usize))
                .collect()
        } else {
            Vec::new()
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:414-426
        for (q, heap) in heaps.iter().enumerate() {
            self.observer
                .heap(q + inactive_heaps.len(), "INIT", 0.0, 0, -1, 0, heap);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3120-3132

        let mut kappa = if self.options.composition != 0 {
            Some(KappaState::new(batch, &full, self.options)?)
        } else {
            None
        };
        super::x_timing::end(XEnd::Init, x_t);

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3543
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1688-1707
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream_mt_utils.c:145-170
        let mut result_batches = Vec::new();
        if !self.options.gapped || self.options.composition == 0 {
            loop {
                let lists = stream.batch_read_observed(self.observer);
                if lists.is_empty() {
                    break;
                }
                result_batches.push(lists);
            }
        }
        if !self.options.gapped {
            for lists in result_batches {
                for mut list in lists {
                    for h in &mut list.hsps {
                        let p = &full[h.hsp.context];
                        h.bit_score = (h.hsp.score as f64 * p.ungapped.lambda
                            - p.ungapped.k.glibc_ln())
                            / std::f64::consts::LN_2;
                    }
                    queries[list.query_index].update_observed(list, self.observer);
                }
            }
        } else if let Some(state) = &mut kappa {
            let mut redo_fetched = None;
            x_t_finish = self.kappa_redo_batch(KappaBatch {
                batch,
                full: &full,
                state,
                stream: &mut stream,
                queries: &mut queries,
                heaps: &mut heaps,
                inactive_heaps: &mut inactive_heaps,
                fetched: (&fetched_reads, &fetched_order, &mut redo_fetched),
            })?;
            fetched_utf8 = redo_fetched;
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1640-1707
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1500-1508
            for lists in result_batches {
                let encoded = self.subjects.read_searched(lists[0].oid, self.options)?;
                for list in lists {
                    let input: Vec<_> = list.hsps.iter().map(|h| h.hsp.clone()).collect();
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-407
                    let redone = super::traceback::ordinary_traceback_observed(
                        &input,
                        batch,
                        &full,
                        self.options,
                        &encoded,
                        &mut |h, contained| self.observer.containment(h, contained),
                    )?;
                    let (output, _) = finalize_redone(
                        &redone,
                        batch,
                        &full,
                        self.options,
                        &encoded,
                        self.db_length,
                        list.query_index,
                        list.oid,
                        self.observer,
                    )?;
                    if !output.hsps.is_empty() {
                        queries[list.query_index].update_observed(output, self.observer);
                    }
                }
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:849-865
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1717-1723
        if self.options.culling_limit > 0 {
            super::culling::traceback_pipe(
                &mut queries,
                batch.contexts.iter().map(|c| c.length as i32).collect(),
                self.options.culling_limit,
                prelim_hitlist_size(self.options),
                self.observer,
            );
        }
        for query in &mut queries {
            for list in &mut query.lists {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2049-2070
                if self.options.max_hsps > 0 {
                    self.observer.hsps("MAX_HSPS_IN", list.oid, &list.hsps);
                    trim_by_max_hsps(&mut list.hsps, self.options.max_hsps as usize);
                    self.observer.hsps("MAX_HSPS_OUT", list.oid, &list.hsps);
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2015-2045
                // The context of an HSP is its query frame, whose length is the frame's.
                if self.options.query_cov_hsp_perc > 0.0
                    && reap_by_query_coverage(
                        &mut list.hsps,
                        self.options.query_cov_hsp_perc,
                        |h| query_extent(&h.hsp, &batch.contexts),
                    )
                {
                    list.refresh();
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2538-2588
                // NCBI calls Blast_HSPListSubjectBestHit on a list the coverage left empty after
                // freeing it (blast_traceback.c:853-861); LOSAT skips the empty list.
                if self.options.subject_besthit && !list.hsps.is_empty() {
                    self.observer.hsps("SUBJECT_BEST_IN", list.oid, &list.hsps);
                    subject_besthit(&mut list.hsps);
                    self.observer.hsps("SUBJECT_BEST_OUT", list.oid, &list.hsps);
                }
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:864-866
        if self.options.query_cov_hsp_perc > 0.0 {
            for query in &mut queries {
                query.lists.retain(|list| !list.hsps.is_empty());
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1777
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1777
        // Not where `BlastSeqSrcGetTotLen` (`Subjects::total_length`) is 0 or negative: the lists
        // keep the traceback's order, by OID through `update_observed` (and the composition
        // heap's reversed order with composition-based statistics), NCBI's one-thread order.
        if self.subjects.total_length()? > 0 {
            for query in &mut queries {
                query.lists.sort_by(super::results::compare_list);
            }
        }
        for (q, query) in queries.iter().enumerate() {
            self.observer.targets("IN", q, &query.lists);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:878-894
        for query in &mut queries {
            query.lists.truncate(self.options.hitlist_size as usize);
        }
        for (q, query) in queries.iter().enumerate() {
            self.observer.targets("OUT", q, &query.lists);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-101
        self.output.push(BatchResults {
            query_ordinal: self.ordinal,
            queries: queries.into_iter().map(|q| q.lists).collect(),
            prepared: batch.x_for_results(),
            parameters: full,
            search_skipped: false,
            fetched_utf8,
            request_space: None,
            unlisted: Vec::new(),
            shown: Vec::new(),
            gapped_k: None,
        });
        super::x_timing::end(XEnd::Finish, x_t_finish);
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3670-3703
#[allow(clippy::too_many_arguments)]
fn finalize_redone(
    redone: &[super::kappa::RedoneHsp],
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    db_length: i64,
    query_index: usize,
    oid: usize,
    observer: &mut dyn RuntimeObserver,
) -> Result<(HspList, i32)> {
    let input: Vec<_> = redone.iter().map(|h| h.hsp.clone()).collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
    if !input.is_empty() {
        if let Some(link) = super::linking::link_parameters(options) {
            observer.link_run(subject.len() as i32, true, &link);
        }
    }
    let mut linked = evaluate_gapped(
        &input,
        batch,
        params,
        options,
        subject.len() as i32,
        db_length,
        if options.composition != 0 { 32.0 } else { 1.0 },
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1976-2008
    let observe_reap = !linked.hsps.is_empty();
    if observe_reap {
        observer.reap("FINAL_REAP_IN", oid, &linked, options.evalue);
    }
    reap_final(&mut linked, options);
    if observe_reap {
        observer.reap("FINAL_REAP_OUT", oid, &linked, options.evalue);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:438-444
    let best_score = linked.hsps.first().map(|h| h.hsp.score).unwrap_or(0);
    let bits = normalize_scores(&mut linked, options);
    let hsps = linked
        .hsps
        .iter()
        .zip(bits)
        .map(|(h, bit_score)| {
            let r = &redone[h.source_index];
            let (identity, positive) = identity_and_positive(r, batch, subject, options.scoring);
            Hsp {
                hsp: h.hsp.clone(),
                num: h.num,
                evalue: h.evalue,
                bit_score,
                identity,
                positive,
                edit_script: r.edit_script.clone(),
                composition_method: r.composition_method,
            }
        })
        .collect();
    let list = HspList {
        query_index,
        oid,
        hsps,
        best_evalue: linked.best_evalue,
    };
    Ok((list, best_score))
}
