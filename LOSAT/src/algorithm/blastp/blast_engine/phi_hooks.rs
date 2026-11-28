//! psiblast `-phi_pattern` (item PS8a): round 1 as a PHI-BLAST search in blastp's call
//! (`run_call`), with the stages of `algorithm::phi_search` where blastp's differ. The pattern
//! search takes the place of the word lookup and the ungapped stage: the queries' contexts without
//! a word table (`phi_lookup`), each query's setup (`phi_start`: the Karlin block of NCBI's table,
//! the occurrences, the search space and the cutoff), the preliminary stage over the subjects
//! (`phi_prelim`: every occurrence against every subject hit, kept by blastp's preliminary hit
//! lists), and the traceback of the kept lists once every subject is scanned (`phi_traceback`: the
//! e-values read the subjects' count of pattern hits). The final lists are blastp's; after
//! `finish_hit_lists` and `report_values`, `PhiCall::split_occurrences` splits them by query
//! occurrence (`PHIBlast_HSPResultsSplit`, `split_hit_list_by_occurrence`) using the typed
//! occurrence index (`PhiOccurrenceIndex`) carried on `BlastpHsp` without growing `BlastpHsp`
//! (168 B). Every function is out of line and cold: blastp's default path only tests
//! `ResolvedBlastpArgs::phi`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1775-1778

use anyhow::Result;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::Mutex;

use super::{
    database, BlastpRuntimeLookup, BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch,
    PrelimCollector, PrelimFilter,
};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{
    update_best_evalue, BlastpHitList, BlastpHsp, BlastpHspList, PhiOccurrenceIndex,
};
use crate::algorithm::blastp::x_prelim_compact::{PhiPrelimOrdinal, XPrelimHsp};
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::phi_search::engine::{self, PhiHsp, PhiSearch, PhiSetup, PhiSetupError};
use crate::algorithm::phi_search::pattern::PhiError;
use crate::algorithm::phi_search::stats;
use crate::algorithm::pssm_search::gapalign::PssmGapScratch;
use crate::algorithm::pssm_search::lookup::build_pssm_lookup_in;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::algorithm::tblastx::lookup::{
    build_lookup_segments_for_context, prepare_matrix_lookup_query_for_word_size, QueryContext,
};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::blastinput::seq_range;

thread_local! {
    static PHI_SCRATCH: RefCell<PssmGapScratch> = RefCell::new(PssmGapScratch::new());
}

/// The queries' contexts of a PHI-BLAST search (the letters that the search reads, as blastp's
/// query preparation makes them) and a word table without words: the pattern scan finds the seeds
/// (`phi_prelim`), and the table is what the call's other readers of the scan's shape take. Every
/// context is valid: PHI-BLAST's score block is the pattern table's (`s_PHIScoreBlkFill`), so
/// NCBI computes no context's Karlin-Altschul parameters from its letters, and a query that blastp
/// cannot score (all `X`) is searched; a context without letters is not valid
/// (`call_setup::call_lookup`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:121-130
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:504-506
#[cold]
#[inline(never)]
pub(super) fn phi_lookup(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    let (_, _, mut contexts) = prepare_matrix_lookup_query_for_word_size(
        query_frames,
        args.word_size,
        false,
        args.scoring.matrix,
        None,
    );
    for context in &mut contexts {
        context.is_valid = true;
    }
    let table = build_pssm_lookup_in(&QueryPssm { rows: Vec::new() }, 3, 0, &[]);
    (BlastpRuntimeLookup::Aa(table), contexts)
}

/// The Karlin block of NCBI's PHI-BLAST table for the search's matrix and gap costs
/// (`stats::phi_karlin`), or NCBI's message where the table has none: NCBI then writes the message
/// in the report of each query and searches nothing (`PhiCall::no_statistics`). The table's names
/// are compared with the matrix name as given (`strcmp`): `-matrix blosum62` has no entry (probe
/// `k.lpam.0`: `Matrix pam250 not allowed in PHI-BLAST`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:165-308
#[cold]
#[inline(never)]
pub(super) fn phi_karlin(args: &ResolvedBlastpArgs) -> Result<stats::PhiKarlin, String> {
    stats::phi_karlin(
        &args.matrix_name,
        args.scoring.gap_open,
        args.scoring.gap_extend,
    )
}

/// The searches of the call's queries and what the preliminary stage gathers for the traceback.
pub(super) struct PhiCall {
    /// Each query's search, by query: `None` for a query without letters (its batch is searched
    /// as its query with letters, `PhiBatch`), and for every query where NCBI's table has no
    /// Karlin block (`no_statistics`). A query whose letters give blastp no Karlin block (all `X`)
    /// is searched: PHI-BLAST's score block is the table's (`s_PHIScoreBlkFill`).
    pub(super) searches: Vec<Option<PhiSearch>>,
    /// The call's query batches, each the result set of one report (`PhiBatch`).
    pub(super) batches: Vec<PhiBatch>,
    /// By query, the letters that the reports read for a query whose results NCBI reports under
    /// its batch's first query without letters (`PhiBatch::shows_gaps`): a gap character at each
    /// of its positions; `None` for the others (`ShownQueries::of`).
    pub(super) shown_gaps: Vec<Option<Box<[u8]>>>,
    /// The pattern hits of the scanned subjects (NCBI's `Int8` `lookup_hits`; the traceback reads
    /// them as an `Int4`, `num_patterns_db`).
    lookup_hits: AtomicI64,
    /// LOSAT's refusal of a subject's scan (`phi_search::scan`), the first one met.
    refusal: Mutex<Option<String>>,
    /// NCBI's message where its table has no Karlin block for the matrix and gap costs
    /// (`phi_karlin`): no query is searched (every `searches` entry is `None`), and each batch's
    /// report is the message, no hits and the footer of -1 values (`s_PHIScoreBlkFill` leaves the
    /// score block without Karlin blocks; NCBI's setup goes on, so the pattern's and the
    /// occurrences' errors still stop the run).
    pub(super) no_statistics: Option<String>,
    /// Per-query, per-occurrence final hit lists (`PHIBlast_HSPResultsSplit`), populated by
    /// `split_occurrences` after `shift_hsp_to_records`.
    occurrence_lists: Vec<Vec<Option<BlastpHitList>>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1652-1712
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:263-281
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:303-360
/// A query batch of PHI-BLAST's call: the queries without letters that psiblast's reader puts
/// before the next query with letters, and that query (`read_psiblast_query_batches` with
/// CLocalBlast's rule). NCBI searches the batch as one query (a context without letters shares
/// the offset of the next, so the occurrences are the lettered query's own), puts every HSP in
/// context 0, and makes each occurrence's result with the batch's first query: its Seq-id, its
/// length and Bioseq for the rows (`PhiBlastResults2SeqAlign_OMF`), the statistics of its context
/// for the footer (`BuildBlastAncillaryData`), and result `i` carries the messages of the batch's
/// query `i` (`BlastBuildSearchResultSet`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct PhiBatch {
    /// The batch's queries; the first heads its report.
    pub(super) queries: std::ops::Range<usize>,
    /// The batch's query with letters (its last), whose search and hit lists the report shows.
    pub(super) searched: usize,
}

impl PhiBatch {
    /// The query that heads the batch's report.
    pub(super) fn first(&self) -> usize {
        self.queries.start
    }

    /// Whether the report shows the rows of `searched` with the letters of a first query without
    /// letters: a gap character at each position (`PhiCall::shown_gaps`).
    pub(super) fn shows_gaps(&self) -> bool {
        self.searched != self.queries.start
    }
}

/// The batches of the call's queries `query_batches` (call-local), each with its query with
/// letters (`PhiBatch`): the last query of the batch that has letters, or its last query.
fn phi_batches(
    query_batches: &[seq_range::RangedBatch],
    query_frames: &[Vec<QueryFrame>],
) -> Vec<PhiBatch> {
    query_batches
        .iter()
        .filter(|batch| !batch.searched.is_empty())
        .map(|batch| PhiBatch {
            queries: batch.searched.clone(),
            searched: batch
                .searched
                .clone()
                .rev()
                .find(|&query| query_frames[query][0].aa_len > 0)
                .unwrap_or(batch.searched.end - 1),
        })
        .collect()
}

impl PhiCall {
    /// The batch whose report query `q_idx` heads, `None` for a query that heads none.
    pub(super) fn batch_at(&self, q_idx: usize) -> Option<&PhiBatch> {
        let at = self
            .batches
            .binary_search_by_key(&q_idx, PhiBatch::first)
            .ok()?;
        Some(&self.batches[at])
    }

    /// The batch's search, `None` where nothing is searched (`no_statistics`).
    pub(super) fn search_of(&self, batch: &PhiBatch) -> Option<&PhiSearch> {
        self.searches.get(batch.searched).and_then(Option::as_ref)
    }

    /// The query of each result of the batch, whose Seq-id and title the result's report shows:
    /// the batch's first for each occurrence of its search (`PhiBlastResults2SeqAlign_OMF`), or,
    /// where nothing is searched (`no_statistics`), each query of the batch for a result of its own
    /// without a Seq-align set (NCBI BLAST+ 2.17.0 `PS/tmp/fix8a/ncbi5/er.first.p250.7`).
    pub(super) fn result_queries(&self, batch: &PhiBatch) -> Vec<usize> {
        match self.search_of(batch) {
            Some(search) => vec![batch.first(); search.occurrences.len()],
            None => batch.queries.clone().collect(),
        }
    }

    /// The footer of the batch's report (`x_PrintOneQueryFooter` of its first query's context).
    pub(super) fn footer(
        &self,
        batch: &PhiBatch,
    ) -> crate::report::pairwise::psiblast_round::phi_round::PhiFooter {
        use crate::report::pairwise::psiblast_round::phi_round::PhiFooter;
        if self.no_statistics.is_some() {
            PhiFooter::NoStatistics
        } else if batch.shows_gaps() {
            PhiFooter::NoLetters
        } else {
            PhiFooter::Searched
        }
    }

    /// The messages that the batch's report writes before it (`PrintPhiResult` writes each
    /// result's warnings): result `i` carries the messages of the batch's query `i`, written as
    /// `per_query` holds them (`query_warnings::phi_result_warning`).
    pub(super) fn result_warnings(&self, batch: &PhiBatch, per_query: &[Vec<u8>]) -> Vec<u8> {
        let results = self.result_queries(batch).len();
        per_query[batch.queries.clone()]
            .iter()
            .take(results)
            .flatten()
            .copied()
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3572-3632
    /// Splits every searched query's final hit list into one hit list per pattern occurrence
    /// (`PHIBlast_HSPResultsSplit`), consuming each HSP's `PhiOccurrenceIndex` attached at
    /// `phi_traceback` (and resetting `subject_frame` to `0`).
    /// Called after `finish_hit_lists`, `report_values`, and `shift_hsp_to_records`.
    #[cold]
    #[inline(never)]
    pub(super) fn split_occurrences(&mut self, hit_lists: &mut [Option<BlastpHitList>]) {
        self.occurrence_lists = (0..self.searches.len())
            .map(|q_idx| {
                let num_occurrences = self.searches[q_idx]
                    .as_ref()
                    .map_or(0, |search| search.occurrences.len());
                split_hit_list_by_occurrence(
                    hit_lists.get_mut(q_idx).and_then(Option::as_mut),
                    num_occurrences,
                    |_oid, hsp| {
                        hsp.take_phi_occurrence()
                            .expect("every PHI HSP carries a PhiOccurrenceIndex from phi_traceback")
                            .index()
                    },
                )
            })
            .collect();
    }

    /// The per-occurrence final hit lists of query `q_idx` (`split_occurrences`).
    pub(super) fn occurrence_lists(&self, q_idx: usize) -> &[Option<BlastpHitList>] {
        self.occurrence_lists.get(q_idx).map_or(&[], Vec::as_slice)
    }
}

/// How the call's PHI-BLAST setup ends.
pub(super) enum PhiStart {
    /// Every query is set up.
    Search(PhiCall),
    /// NCBI's error of a query's setup: the pattern's compilation (`query` 0, the first query that
    /// NCBI sets up) or the query's occurrences; NCBI's engine stops with it after the reports of
    /// the queries before (`BLAST engine error: Warning: <message> `).
    EngineError { query: usize, text: String },
}

/// What LOSAT's limits of PHI-BLAST read besides the options (`phi_limits`): the call's formats
/// and its subjects.
pub(super) struct PhiLimits<'a> {
    formats: &'a [super::BlastpReport],
    subject_source: super::BlastpSubjectSource<'a>,
}

impl<'a> PhiLimits<'a> {
    /// The limits' inputs: the call's formats and subjects.
    pub(super) fn new(
        formats: &'a [super::BlastpReport],
        subject_source: super::BlastpSubjectSource<'a>,
    ) -> Self {
        Self {
            formats,
            subject_source,
        }
    }
}

/// psiblast `-phi_pattern`: the run's PHI-BLAST call (`phi_start` with the table's Karlin block,
/// `phi_karlin`), `None` without a pattern; NCBI's engine error of a query's setup is
/// `setup_error(query, text)` (`phi_setup_error` of `run_call`). LOSAT's limits (`phi_limits`)
/// come after NCBI's errors of the first batch's setup (the pattern's, its query's occurrences)
/// and after the table's lookup, where NCBI would search what LOSAT does not: before an error of a
/// later batch, whose batches before NCBI searches and reports.
pub(super) fn phi_call(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    query_batches: &[seq_range::RangedBatch],
    limits: PhiLimits<'_>,
    db_length: i64,
    db_num_seqs: i64,
    setup_error: impl FnOnce(usize, &str) -> anyhow::Error,
) -> Result<Option<PhiCall>> {
    if args.phi.is_none() {
        return Ok(None);
    }
    let karlin = phi_karlin(args);
    let searched = karlin.is_ok();
    let batches = phi_batches(query_batches, query_frames);
    let start = phi_start(args, query_frames, batches, karlin, db_length, db_num_seqs)?;
    let first_batch_fails = match &start {
        PhiStart::EngineError { query, .. } => query_batches
            .iter()
            .find(|batch| !batch.searched.is_empty())
            .is_none_or(|first| first.searched.contains(query)),
        PhiStart::Search(_) => false,
    };
    if !first_batch_fails {
        if let Some(refusal) = phi_limits(args, &limits, searched) {
            return Err(refusal);
        }
    }
    match start {
        PhiStart::Search(call) => Ok(Some(call)),
        PhiStart::EngineError { query, text } => Err(setup_error(query, &text)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-900
/// LOSAT's limits of PHI-BLAST (item PS8a; their order after NCBI's errors is `phi_call`'s):
/// outfmt 5 and 11 of the round (a later item, PS8b) and `-html` (NCBI's HTML PHI report is not
/// ported yet), whatever is searched. Then only where the
/// round is searched (`searched`: NCBI's table has a Karlin block for the matrix and gap costs;
/// without it round 1 is the run's last whatever `-num_iterations` and `-sum_stats` say, probes
/// `k.pam.it3.0`, `k.pam.ss3.0`): the rounds after it (PS8b), `-sum_stats true` (NCBI links the
/// preliminary HSPs with the pattern's Karlin block and reaps them by the linked e-values; probe
/// `op.sumstats.6`) and a subject longer than `MAX_DBSEQ_LEN` (NCBI's parts), both item PS8a-2.
/// Outfmt 8, 9 and 12 are no exception: NCBI's engine error of their Seq-annot comes only where an
/// alignment is found (`report_write::write_phi_report`; `-sum_stats true` with outfmt 8 ends with
/// exit status 0 and no output, probe `PS/tmp/fix8a/ncbi5/k.ss.8`).
fn phi_limits(
    args: &ResolvedBlastpArgs,
    limits: &PhiLimits<'_>,
    searched: bool,
) -> Option<anyhow::Error> {
    let program = args.app.name();
    use super::BlastpReport;
    if let Some(number) = limits.formats.iter().find_map(|format| match format {
        BlastpReport::Xml => Some(5),
        BlastpReport::Archive => Some(11),
        _ => None,
    }) {
        return Some(anyhow::anyhow!(
            "-phi_pattern with output format {number} is not supported by LOSAT's {program} yet"
        ));
    }
    if args.html.is_some() {
        return Some(anyhow::anyhow!(
            "-phi_pattern with -html is not supported by LOSAT's {program} yet"
        ));
    }
    if !searched {
        return None;
    }
    if let Some(rounds) = args.psi_rounds.as_deref() {
        return Some(anyhow::anyhow!(
            "-phi_pattern with -num_iterations {} is not supported by LOSAT's {program} yet",
            rounds.num_iterations
        ));
    }
    if args.sum_stats {
        return Some(anyhow::anyhow!(
            "-phi_pattern with -sum_stats true is not supported by LOSAT's {program} yet"
        ));
    }
    let (what, length) = limits.subject_source.longest_split()?;
    Some(super::too_long_phi_subject_refusal(&what, length))
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:507-532
/// What psiblast sets in the options handle that the search strategy's requests hold
/// (`report::archive::OptionsHandle::set_search_sets`): after each batch it resets the
/// composition-based statistics to their value before the batches, which PHI-BLAST's options never
/// set (psiblast's arguments do not set them for `phiblastp`), so the reset appends
/// `CompositionBasedStats 0` after the batch's `EffectiveSearchSpace`: the requests of the batches
/// after the first hold both (NCBI BLAST+ 2.17.0 `st.multi`).
#[cold]
#[inline(never)]
pub(super) fn note_handle_sets(handle: &crate::report::archive::OptionsHandle) -> Result<()> {
    use crate::report::archive::fields;
    use crate::report::archive::options::{Param, ParamValue};
    handle.set_search_sets(
        Vec::new(),
        Vec::new(),
        vec![Param::new(
            fields::COMPOSITION_BASED_STATS,
            ParamValue::Integer(0),
        )],
    )
}

/// The setup of each query of the call: NCBI's errors of the pattern (`PhiError::Ncbi`) and of
/// the occurrences of each batch's query with letters (`PhiSetupError`), LOSAT's refusals of the
/// query's scan, or the searches. A query without letters is not set up (its batch's query with
/// letters is, `PhiBatch`). `db_length` and `db_num_seqs` are the database's (`-dbsize` replaces
/// the length). The Karlin block of NCBI's table, or its message where there is none, is the
/// caller's (`phi_karlin`); without it the queries' occurrences are found for their errors alone.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1065-1118
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:170-185
#[cold]
#[inline(never)]
pub(super) fn phi_start(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    batches: Vec<PhiBatch>,
    karlin: Result<stats::PhiKarlin, String>,
    db_length: i64,
    db_num_seqs: i64,
) -> Result<PhiStart> {
    let phi = args
        .phi
        .as_deref()
        .expect("the setup of a PHI-BLAST search");
    let pattern = match phi {
        Ok(pattern) => pattern,
        Err(PhiError::Ncbi(message)) => {
            return Ok(PhiStart::EngineError {
                query: 0,
                text: crate::report::query_warnings::setup_error(message),
            });
        }
        // `validate_requested_blastp_support` has refused it.
        Err(PhiError::Refusal(message)) => anyhow::bail!("{message}"),
    };
    let rows = crate::utils::matrix::protein_score_rows(args.scoring.matrix);
    let score = |q: u8, s: u8| rows.score(q, s);
    let mut searches = Vec::with_capacity(query_frames.len());
    for (query, frames) in query_frames.iter().enumerate() {
        let frame = &frames[0];
        // A query without letters: its batch is its query with letters' search.
        if frame.aa_len == 0 {
            searches.push(None);
            continue;
        }
        let letters = &frame.aa_seq[1..=frame.aa_len];
        let segments: Vec<(usize, usize)> =
            build_lookup_segments_for_context(0, frame.aa_len, &frame.seg_masks)
                .into_iter()
                .map(|(left, right)| (left as usize, right as usize))
                .collect();
        // Without NCBI's Karlin block the query is set up for its errors alone.
        let set_up = match &karlin {
            Ok(karlin) => engine::setup(&PhiSetup {
                pattern: std::sync::Arc::clone(pattern),
                query: letters,
                segments: &segments,
                score: &score,
                karlin: *karlin,
                db_length,
                db_num_seqs,
                evalue: args.evalue,
                gap_open: args.scoring.gap_open,
                gap_extend: args.scoring.gap_extend,
                x_drop_bits: args.x_drop_gapped_bits,
                x_drop_final_bits: args.x_drop_final_bits,
            })
            .map(Some),
            Err(_) => engine::occurrences(pattern, letters, &segments).map(|_| None),
        };
        match set_up {
            Ok(search) => searches.push(search),
            Err(PhiSetupError::Refusal(message)) => anyhow::bail!("{message}"),
            Err(error) => {
                return Ok(PhiStart::EngineError {
                    query,
                    text: crate::report::query_warnings::setup_error(&String::from_utf8_lossy(
                        &error.message(),
                    )),
                });
            }
        }
    }
    let mut shown_gaps = vec![None; query_frames.len()];
    for batch in batches.iter().filter(|batch| batch.shows_gaps()) {
        shown_gaps[batch.searched] =
            Some(vec![b'-'; query_frames[batch.searched][0].aa_len].into_boxed_slice());
    }
    Ok(PhiStart::Search(PhiCall {
        searches,
        batches,
        shown_gaps,
        lookup_hits: AtomicI64::new(0),
        refusal: Mutex::new(None),
        no_statistics: karlin.err(),
        occurrence_lists: Vec::new(),
    }))
}

/// The preliminary HSP that the hit lists keep for a PHI HSP of query `q_idx` and subject
/// `s_idx`: its extent and score (the e-value of the preliminary stage is 0, NCBI computes none).
fn preliminary_of(hsp: &PhiHsp, q_idx: usize, s_idx: usize, query_length: usize) -> XPrelimHsp {
    XPrelimHsp {
        hsp: BlastpPreliminaryHsp {
            query_start: hsp.query_start,
            query_end: hsp.query_end,
            subject_start: hsp.subject_start,
            subject_end: hsp.subject_end,
            gapped_query_start: hsp.q_pattern_offset,
            gapped_subject_start: hsp.s_pattern_offset,
            raw_score: hsp.score,
            query_context: q_idx as i32,
            query_frame: 0,
            subject_frame: 0,
            query_length,
            q_idx: q_idx as u32,
            s_idx: s_idx as u32,
        },
        e_value: 0.0,
    }
}

/// The preliminary stage of a PHI-BLAST search (`PHIBlastWordFinder` and `PHIGetGappedScore`):
/// each subject's pattern hits against every query occurrence (`PhiSearch::prelim`), the HSPs at
/// the cutoff or above kept by the preliminary hit lists (blastp's collector, or the writer of
/// `-culling_limit` / `-best_hit_*`), and the subject's hits counted for the traceback. A
/// database is read in blocks as blastp's (`database::search_blocks`), a subject set in turn.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1016-1046
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn phi_prelim(
    phi: &PhiCall,
    contexts: &[QueryContext],
    window: i32,
    new_scratch: &(dyn Fn() -> BlastpSubjectScratch + Sync),
    database: Option<&database::BlastpDatabase>,
    subjects: &[EncodedProtein],
    pool: &crate::utils::threading::SearchPool<'_>,
    prelim_filter: Option<&mut PrelimFilter>,
    prelim_hitlist_size: usize,
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) -> Result<()> {
    // Without NCBI's Karlin block nothing is searched (`PhiCall::no_statistics`).
    if phi.searches.iter().all(Option::is_none) {
        return Ok(());
    }
    let process_subject = |_: &mut BlastpSubjectScratch,
                           s_idx: usize,
                           subject: BlastpSubjectRef<'_>,
                           subject_result: &mut BlastpSubjectResult| {
        let letters = &subject.aa_seq[1..];
        let mut counted = false;
        for (q_idx, search) in phi.searches.iter().enumerate() {
            // A query without letters: its batch's query with letters is searched (`PhiBatch`).
            let Some(search) = search else {
                continue;
            };
            let found = PHI_SCRATCH
                .with(|scratch| search.prelim(letters, subject.aa_len, &mut scratch.borrow_mut()));
            let (hsps, hits) = match found {
                Ok(found) => found,
                Err(PhiError::Refusal(message)) | Err(PhiError::Ncbi(message)) => {
                    let mut refusal = phi.refusal.lock().unwrap_or_else(|e| e.into_inner());
                    refusal.get_or_insert(message);
                    return;
                }
            };
            // The hits are the subject's, whatever the query: counted once per subject.
            if !counted {
                phi.lookup_hits.fetch_add(hits as i64, Ordering::Relaxed);
                counted = true;
            }
            if hsps.is_empty() {
                continue;
            }
            let query_length = contexts[q_idx].aa_len;
            let list = BlastpHspList {
                oid: s_idx as u32,
                query_index: q_idx as u32,
                hsps: hsps
                    .iter()
                    .enumerate()
                    .map(|(index, hsp)| {
                        preliminary_of(hsp, q_idx, s_idx, query_length)
                            .with_phi_ordinal(PhiPrelimOrdinal::new(index))
                    })
                    .collect(),
                best_evalue: 0.0,
            };
            subject_result.x_compact_hsp_lists.push((q_idx, list));
        }
    };
    let searched = match (database, prelim_filter) {
        (Some(database), None) => database::search_blocks(
            database,
            pool,
            window,
            &new_scratch,
            &process_subject,
            &mut PrelimCollector {
                x_compact_hit_lists,
                prelim_hitlist_size,
            },
        ),
        (Some(database), Some(filter)) => database::search_blocks(
            database,
            pool,
            window,
            &new_scratch,
            &process_subject,
            filter,
        ),
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1498
        (None, mut prelim_filter) => {
            let mut scratch = new_scratch();
            let mut subject_result = BlastpSubjectResult::new(0);
            for (s_idx, subject) in subjects.iter().enumerate() {
                process_subject(
                    &mut scratch,
                    s_idx,
                    BlastpSubjectRef::of(subject),
                    &mut subject_result,
                );
                if let Some(filter) = prelim_filter.as_deref_mut() {
                    filter.add(&mut subject_result);
                    continue;
                }
                for (q_idx, list) in subject_result.x_compact_hsp_lists.drain(..) {
                    x_compact_hit_lists[q_idx]
                        .get_or_insert_with(|| BlastpHitList::new(prelim_hitlist_size))
                        .update(list);
                }
            }
            Ok(())
        }
    };
    searched?;
    if let Some(message) = phi.refusal.lock().unwrap_or_else(|e| e.into_inner()).take() {
        anyhow::bail!("{message}");
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:780-794
/// Keeps only the re-scanned `PhiHsp`s that survived in the preliminary hit list `kept`
/// (after `-culling_limit` / `-best_hit_*` preliminary filtering), selecting each entry by its
/// `PhiPrelimOrdinal` in `kept` order so variable-length pattern hits sharing a start offset and
/// preliminary bounds stay unambiguous.
fn select_kept_prelim(kept: &[XPrelimHsp], all: Vec<PhiHsp>) -> Vec<PhiHsp> {
    let mut slots: Vec<Option<PhiHsp>> = all.into_iter().map(Some).collect();
    kept.iter()
        .map(|prelim| {
            let index = prelim
                .phi_ordinal()
                .expect("every PHI preliminary HSP carries a PhiPrelimOrdinal")
                .index();
            slots
                .get_mut(index)
                .and_then(Option::take)
                .expect("each PHI preliminary ordinal is in range and unique")
        })
        .collect()
}

/// The traceback of a PHI-BLAST search (`s_PHITracebackFromHSPList` of each kept list): the
/// subject's preliminary HSPs that survived in `list.hsps` (re-found by `PhiSearch::prelim` and
/// selected by `select_kept_prelim` via `PhiPrelimOrdinal`), each traced with the final X-drop,
/// its e-value from the count of pattern hits of every scanned subject (`num_patterns_db`, NCBI's
/// `(Int4)` of its `Int8` sum), reaped by `-evalue`; the final HSPs with their counts
/// (`pssm_search::traceback0::final_hsp`), into the final hit lists as blastp's traceback inserts
/// them (`traceback::hit_lists_of_traced`). Each final HSP carries its `PhiOccurrenceIndex`
/// (`pat_info->index`, `Int4` range) for `PhiCall::split_occurrences`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:780-815
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1636-1641
#[cold]
#[inline(never)]
pub(super) fn phi_traceback(
    phi: &mut PhiCall,
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    subjects: XRedoSubjects<'_>,
    compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    sorts_results: bool,
) -> Result<Vec<Option<BlastpHitList>>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1775-1778
    let num_patterns_db = phi.lookup_hits.load(Ordering::Relaxed) as i32;
    let mut traced: Vec<(usize, Option<BlastpHspList>)> = Vec::new();
    let mut scratch = PssmGapScratch::new();
    for (q_idx, hit_list) in compact_hit_lists.iter_mut().enumerate() {
        let Some(hit_list) = hit_list.take() else {
            continue;
        };
        let search = phi.searches[q_idx]
            .as_ref()
            .expect("the hit list of a searched query");
        let ctx = &contexts[q_idx];
        let query_nomask = super::query_nomask_sequence(ctx);
        let count = hit_list.hsplist_count;
        for list in hit_list.hsplist_array.into_iter().take(count) {
            if list.hsps.is_empty() {
                traced.push((q_idx, None));
                continue;
            }
            let read = subjects.read(list.oid)?;
            let residues = read.residues();
            let mut subject = Vec::with_capacity(residues.len() + 1);
            subject.extend_from_slice(residues);
            subject.push(0);
            let (all_hsps, _) = search
                .prelim(&subject, residues.len(), &mut scratch)
                .map_err(|error| match error {
                    PhiError::Refusal(message) | PhiError::Ncbi(message) => {
                        anyhow::anyhow!("{message}")
                    }
                })?;
            let hsps = select_kept_prelim(&list.hsps, all_hsps);
            let kept = search
                .traceback(
                    &subject,
                    residues.len(),
                    hsps,
                    num_patterns_db,
                    args.evalue,
                    &mut scratch,
                )
                .map_err(|error| match error {
                    PhiError::Refusal(message) | PhiError::Ncbi(message) => {
                        anyhow::anyhow!("{message}")
                    }
                })?;
            if kept.is_empty() {
                traced.push((q_idx, None));
                continue;
            }
            let hsps = kept
                .into_iter()
                .map(|hsp| {
                    let mut preliminary =
                        preliminary_of(&hsp, q_idx, list.oid as usize, ctx.aa_len).hsp;
                    preliminary.query_end = hsp.query_end;
                    preliminary.subject_end = hsp.subject_end;
                    let mut final_hsp = crate::algorithm::pssm_search::traceback0::final_hsp(
                        &preliminary,
                        hsp.edit_script,
                        query_nomask,
                        residues,
                        args.scoring.matrix,
                    );
                    final_hsp.e_value = hsp.evalue;
                    final_hsp.bit_score = hsp.bit_score;
                    final_hsp.with_phi_occurrence(PhiOccurrenceIndex::new(hsp.occurrence))
                })
                .collect();
            let mut traced_list = BlastpHspList {
                oid: list.oid,
                query_index: list.query_index,
                hsps,
                best_evalue: 0.0,
            };
            update_best_evalue(&mut traced_list);
            traced.push((q_idx, Some(traced_list)));
        }
    }
    Ok(super::traceback::hit_lists_of_traced(
        traced,
        compact_hit_lists.len(),
        args.traceback_hitlist_size(),
        sorts_results,
    ))
}

/// Splits a query's final hit list into one hit list per pattern occurrence
/// (`PHIBlast_HSPResultsSplit`): each surviving HSP is placed in the hit list of its occurrence
/// (`occurrence_of(oid, hsp)`), preserving each subject's HSP order, and each non-empty
/// occurrence hit list is sorted by `s_EvalueCompareHSPLists` (best e-value 0, then top score,
/// then the higher subject first).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3572-3632
#[cold]
#[inline(never)]
pub(super) fn split_hit_list_by_occurrence(
    hit_list: Option<&mut BlastpHitList>,
    num_occurrences: usize,
    mut occurrence_of: impl FnMut(u32, &mut BlastpHsp) -> usize,
) -> Vec<Option<BlastpHitList>> {
    let mut parts: Vec<Option<BlastpHitList>> = vec![None; num_occurrences];
    let Some(hit_list) = hit_list else {
        return parts;
    };
    for list in hit_list
        .hsplist_array
        .iter_mut()
        .take(hit_list.hsplist_count)
    {
        let mut by_occ: BTreeMap<usize, Vec<BlastpHsp>> = BTreeMap::new();
        for hsp in &mut list.hsps {
            let occ = occurrence_of(list.oid, hsp);
            debug_assert!(
                occ < num_occurrences,
                "PHI occurrence {occ} out of range {num_occurrences}"
            );
            if occ < num_occurrences {
                by_occ.entry(occ).or_default().push(hsp.clone());
            }
        }
        for (occ, hsps) in by_occ {
            // NCBI's lists of the split keep the best e-value of a new list, 0
            // (`Blast_HSPListNew`), which the sort below reads: they sort by their first HSP's
            // score, then by subject.
            let part = parts[occ].get_or_insert_with(|| BlastpHitList::new(1));
            part.hsplist_array.push(BlastpHspList {
                oid: list.oid,
                query_index: list.query_index,
                hsps,
                best_evalue: 0.0,
            });
            part.hsplist_count += 1;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3626-3629
    for part in parts.iter_mut().flatten() {
        part.hsplist_max = part.hsplist_count;
        part.sort_by_evalue();
    }
    parts
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:752-761
/// The batch's context 0 as the search strategy's `EffectiveSearchSpace` option reads it
/// (`CallRequests::set_searched`) for the batch whose first query is `first`: the pattern's search
/// space of its search (`PhiSearch::search_space`, the footers' value, `phi_footers`) and the first
/// query's length (0 for a context that is not valid, as `context0_searched` gives it); `None`
/// where nothing is searched (`no_statistics`).
pub(super) fn request_searched(
    phi: &PhiCall,
    contexts: &[QueryContext],
    query_valid: &[bool],
    first: usize,
) -> Option<crate::report::archive::SearchedContext0> {
    let search = phi.search_of(phi.batch_at(first)?)?;
    Some(crate::report::archive::SearchedContext0 {
        space: search.search_space,
        query_length: if query_valid[first] {
            contexts[first].aa_len as i64
        } else {
            0
        },
    })
}

/// The footers of PHI-BLAST's outfmt 0 report: each searched query's ungapped and gapped Karlin
/// blocks are the pattern's (`s_PHIScoreBlkFill`: Lambda, K and H of the table, in both), its
/// search space the pattern's (`PhiSearch::search_space`, which the search strategy's requests
/// carry too, `CallRequests::set_searched`), and the Gumbel columns, where the report prints them
/// (no `-searchsp`), zeros (PHI-BLAST's Gumbel block is not filled). Which footer a batch's report
/// writes is `PhiCall::footer`'s (without the table's Karlin block, the -1 blocks of a query that
/// is not searched; a first query without letters has no Karlin block).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:504-506
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478
pub(super) fn phi_footers(
    phi: &PhiCall,
    queries: &mut [crate::report::BlastpPairwiseQuery],
    report: &mut crate::report::BlastpPairwiseReport,
) {
    let karlin =
        |search: &PhiSearch, base: crate::stats::KarlinParams| crate::stats::KarlinParams {
            lambda: search.karlin.lambda,
            k: search.karlin.k,
            h: search.karlin.h,
            ..base
        };
    for (query, search) in queries.iter_mut().zip(&phi.searches) {
        if let Some(search) = search {
            query.ungapped_karlin = karlin(search, query.ungapped_karlin);
            query.effective_search_space = search.search_space;
        }
    }
    if let Some(search) = phi.searches.iter().flatten().next() {
        report.gapped_karlin = karlin(search, report.gapped_karlin);
    }
    report.gumbel = report.gumbel.map(|gumbel| crate::stats::BlastGumbelBlk {
        a: 0.0,
        alpha: 0.0,
        sigma: 0.0,
        a_un: 0.0,
        alpha_un: 0.0,
        ..gumbel
    });
}

#[cfg(test)]
#[path = "phi_hooks_tests.rs"]
mod tests;
