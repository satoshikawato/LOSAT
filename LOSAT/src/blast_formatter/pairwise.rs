//! blast_formatter's pairwise report (`-outfmt 0`, the default) and query-anchored reports
//! (`-outfmt` 1 to 4), with or without `-html`, of one archive object: NCBI's `CBlastFormat` with `PrintProlog`, the prefetch, `PrintOneResultSet` of
//! each query and `PrintEpilog`, written by the search programs' report writers
//! (`report::pairwise`, `report::pairwise::psiblast_round`) from the archive's Seq-aligns as
//! `CDisplaySeqalign` reads them (decision F3-19 of line RG, O1: the writers' `PairwiseHit`, made
//! here from each Seq-align's Dense-seg and the letters of its query and subject). The
//! query-anchored reports read the sequences, their labels and masks (`ArchiveSequences`), as
//! the search programs' do.
//!
//! The object's scope is the tabular driver's (`tabular`: the request's queries, the database
//! that the archive names or its `-subject` Bioseqs, a query Bioseq found first as NCBI's scope
//! finds it). The report's statistics are the archive's: the options of the first search object
//! (`BuiltOptions`: the factory defaults of its program with the archive's calls, which the
//! epilog prints), the results' first Karlin-Altschul blocks and search space for every query
//! (`report::pairwise::Ancillary::of_archive`), the subjects as the database holds them now or as
//! the request's Bioseqs. What LOSAT does not show as NCBI does is refused by name before the
//! object's first byte (`Refusal::NotCovered`).

use std::cell::RefCell;
use std::collections::HashMap;
use std::io::{self, Write};
use std::ops::Range;
use std::sync::Arc;

use anyhow::Result;

use super::refusals::{Refusal, PAIRWISE_FORMATS};
use super::tabular::{
    covered, formatter_matrix, prefetch_reads, query_sets_taxon, set_titles, BioseqRecords,
    IdIndex, LetterCache, LetterEnd, ObjectInputs, QueryEntry, ReportSink, ScopeBioseqs,
    SubjectSource, Subjects,
};
use super::PROGRAM;
use crate::algorithm::common::database::{
    check_pairwise_parsed, check_reported_title, check_shown_protein_title, pairwise_database,
    pairwise_database_masked, prefetch_messages, prefetch_warning_of, refuse_described_lines_of,
    refuse_described_records_of, ReportedHeader, TitleReads,
};
use crate::blastinput::app::ReportFormat;
use crate::common::{GapEditOp, Hit};
use crate::config::ScoringMatrix;
use crate::report::archive::read::options::{BuiltOptions, EProgram, Setter};
use crate::report::archive::read::request::{Queries, Request};
use crate::report::archive::read::results::{MaskRange, Results};
use crate::report::archive::read::{seq_align_id, Stop};
use crate::report::db_subject::{HitLoad, HitLoads, ParsedSubject};
use crate::report::defline::SourceTaxname;
use crate::report::ostream::ReportWrite;
use crate::report::pairwise::psiblast_round::{
    write_psiblast_epilog, write_psiblast_round_report, ReportRound,
};
use crate::report::pairwise::{
    blastn_section_hsps, blastp_section_aborts, blastp_section_frees, scan_mode_definition_bytes,
    section_aborts, write_blastn_pairwise_report, write_blastp_pairwise_window, Ancillary,
    BlastnIntro, BlastnPairwiseQuery, BlastnPairwiseReport, BlastpPairwiseQuery,
    BlastpPairwiseReport, PairwiseConfig, PairwiseDatabase, PairwiseEdge, PairwiseHit,
    PairwiseView,
};
use crate::report::query_anchored::aln_vec::{AlnMap, MapBuffers};
use crate::report::query_anchored::{fill_seqid, Denseg, Mask, ReportSequences, SeqKey};
use crate::report::query_warnings::QueryWarnings;
use crate::report::serial::seqalign::{Segments, Strand};
use crate::report::serial::seqalign_read::{read_set, PairDenseg, ReadAlign, ReadStop};
use crate::report::tabular_seqalign::RowLetters;
use crate::seqdb::SeqType;
use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::ParsedSeqId;
use crate::stats::tables::KarlinParams;

/// The refusal of a form that the pairwise format does not write yet.
fn refused(what: &str) -> anyhow::Error {
    Refusal::NotCovered {
        format: PAIRWISE_FORMATS,
        what: what.to_string(),
    }
    .into()
}

/// A stop of a Seq-align's read as the run's error: NCBI's where it stops there, else the
/// refusal of the form.
fn read_stop(stop: ReadStop) -> anyhow::Error {
    match stop {
        ReadStop::Stop(stop) => super::archive::stop_error(stop),
        ReadStop::Unsupported(what) => refused(what),
    }
}

/// What the pairwise report reads of the run besides the object (`ObjectInputs`).
pub(super) struct PairwiseInputs<'a> {
    /// The formatting options of the object (`-num_descriptions`, `-num_alignments`,
    /// `-line_length`, `-sorthits`, `-sorthsps`, `-show_gis`, `-html`).
    pub(super) formatting: &'a super::formatting::Formatting,
    /// The options' matrix name when the epilog is printed (`options.GetMatrixName()`: an
    /// `rpsblast` object's `GetResultSet` sets its search statistics' matrix before it), `None`
    /// for NULL.
    pub(super) epilog_matrix: Option<&'a [u8]>,
    /// The masks of each query as `GetResultSet` assigns them (`Results::masks`).
    pub(super) masks: &'a [Vec<MaskRange<'a>>],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:60-178
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:136-219
// NCBI reference (598d8ae6): c++/src/algo/blast/api/disc_nucl_options.cpp:59-89
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:61-99
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:788-813
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:58-157
/// The values of the options that the epilog and the report's tables print
/// (`CBlastFormat::PrintEpilog`, blast_format.cpp:2249-2292): those of the factory options of
/// the program (`CBlastOptionsFactory::Create`: `SetTraditionalBlastnDefaults` for blastn,
/// megablast's, dc-megablast's, rpsblast's handle, which sets no word threshold, and the protein
/// handles of blastp, psiblast, deltablast and PHI-BLAST), then the archive's calls on them in
/// order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct EpilogValues {
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    word_threshold: i32,
    window_size: i32,
    gapped: bool,
    sum_statistics: bool,
}

impl EpilogValues {
    fn of(options: &BuiltOptions) -> Self {
        let (reward, penalty, gap_open, gap_extend, word_threshold, window_size) =
            match options.program {
                EProgram::Blastn => (2, -3, 5, 2, 0, 0),
                EProgram::Megablast => (1, -2, 0, 0, 0, 0),
                EProgram::DiscMegablast => (2, -3, 5, 2, 0, 40),
                EProgram::RpsBlast => (1, -3, 11, 1, 0, 40),
                _ => (1, -3, 11, 1, 11, 40),
            };
        let mut values = EpilogValues {
            reward,
            penalty,
            gap_open,
            gap_extend,
            word_threshold,
            window_size,
            gapped: true,
            sum_statistics: false,
        };
        for setter in &options.setters {
            match *setter {
                Setter::MatchReward(value) => values.reward = value,
                Setter::MismatchPenalty(value) => values.penalty = value,
                Setter::GapOpeningCost(value) => values.gap_open = value,
                Setter::GapExtensionCost(value) => values.gap_extend = value,
                Setter::WordThreshold(value) => values.word_threshold = value,
                Setter::WindowSize(value) => values.window_size = value,
                Setter::GappedMode(gapped) => values.gapped = gapped,
                Setter::SumStatisticsMode(on) => values.sum_statistics = on,
                _ => {}
            }
        }
        values
    }
}

/// The subject of each run of a query's Seq-aligns whose subject ids match
/// (`CDisplaySeqalign`'s and `CShowBlastDefline`'s grouping, showalign.cpp:1972-1975): its
/// source (`Subjects::sources`), by the report's subject index.
#[derive(Default)]
struct RunSubjects {
    sources: Vec<usize>,
}

/// The values of one source as the writers read them (by source index).
struct SourceShown {
    title: Arc<[u8]>,
    parsed: Option<Arc<ParsedSubject>>,
    source: SourceTaxname,
    length: usize,
    /// A database sequence's OID.
    oid: Option<crate::seqdb::Oid>,
}

/// The object's report (`PrintFormattedOutput` of format 0 after the results are read): the
/// prolog, the prefetch, each query's report, the epilog. Every refusal comes before the object's
/// first byte.
pub(super) fn print_object(
    inputs: &ObjectInputs<'_>,
    pairwise: &PairwiseInputs<'_>,
    request: &Request,
    results: &Results,
    sets: &[Range<usize>],
    alignments: &[crate::report::serial::asn_types::view::V<'_>],
    sink: &mut dyn ReportSink,
) -> Result<()> {
    let options = inputs.options;
    if !covered(options.program) {
        return Err(refused(
            "a translated search or a program other than blastn, blastp, psiblast, deltablast, \
             PHI-BLAST of a protein query and rpsblast",
        ));
    }
    if matches!(request.queries, Queries::SeqLocList(_)) {
        return Err(refused(
            "queries given as locations, whose sequences NCBI BLAST+ reads from the database",
        ));
    }
    if set_titles(&request.queries) {
        return Err(refused("a title on the Bioseq-set of the queries"));
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:276-282
    if inputs.database.is_none() && std::env::var_os("BL2SEQ_LEGACY").is_some() {
        return Err(refused(
            "-subject sequences and BL2SEQ_LEGACY set (NCBI BLAST+'s legacy bl2seq report, \
             which LOSAT's search programs refuse too)",
        ));
    }
    let program_type = super::search_program::program_type(options.program);
    let nucleotide = program_type.blastn;
    let protein_subjects = !program_type.subject_nucleotide;
    let program = super::search_program::formatter_program(options.program);
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:452-460
    let psi_task = options.program == EProgram::PsiBlast;
    let values = EpilogValues::of(options);
    if nucleotide && values.word_threshold != 0 {
        return Err(refused(
            "a blastn search with a word threshold, whose epilog line LOSAT's blastn report \
             does not write",
        ));
    }
    if psi_task && matches!(inputs.iteration, Ok(None)) {
        return Err(refused(
            "a psiblast round number of 4294967295, which NCBI BLAST+ reads as no round",
        ));
    }
    let formatting = pairwise.formatting;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:186
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:186-199
    let matrix = formatter_matrix(inputs.matrix_name).unwrap_or(ScoringMatrix::Blosum62);
    let ungapped = !values.gapped;

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:430, 455
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:728-734
    // The prefetch reads `PRE_FETCH_SEQS_LIMIT` after `PrintProlog`: `INT_MAX` loads every
    // subject through NCBI's data loaders (refused before any byte, as for the other formats),
    // a value that is no integer ends the run after the prolog.
    let prefetch_limit = crate::blastinput::app::prefetch_seqs_limit();
    if matches!(prefetch_limit, Ok(Some(i32::MAX))) {
        return Err(Refusal::Prefetch.into());
    }

    let queries: Vec<QueryEntry<'_>> = request
        .queries
        .bioseqs()
        .into_iter()
        .map(|bioseq| QueryEntry::new(bioseq, inputs.believe))
        .collect();
    let subject_bioseqs = request.subject.bioseqs();
    if sets.len() > queries.len() {
        return Err(refused("a result set without its query Bioseq"));
    }
    let scope = ScopeBioseqs {
        query_ids: IdIndex::of_bioseqs(queries.iter().map(|query| query.bioseq)),
        queries: &queries,
        subject_ids: IdIndex::of_bioseqs(subject_bioseqs.iter().copied()),
        subjects: &subject_bioseqs,
    };
    let mut subjects = Subjects::default();

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:658-690, 725-786
    // `BlastFormatter_PreFetchSequenceData` as the tabular driver reads it (`tabular`).
    let seq_aligns: usize = sets.iter().map(|set| set.len()).sum();
    let prefetches = match &prefetch_limit {
        Ok(None) | Ok(Some(i32::MAX)) => true,
        Ok(Some(0)) => false,
        Ok(Some(limit)) => i64::try_from(seq_aligns).unwrap_or(i64::MAX) <= i64::from(*limit),
        Err(_) => false,
    };
    let prefetched = prefetches && prefetch_reads(sets, alignments)?;
    if prefetched {
        for set in sets {
            for &align in &alignments[set.clone()] {
                let id = seq_align_id(align, 1).map_err(super::archive::stop_error)?;
                subjects.index_of(inputs, &scope, id)?;
            }
        }
    }
    let prefetched_subjects = subjects.sources.len();

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:617-656
    let view = ReportFormat::from_number(inputs.choice.number)
        .and_then(PairwiseView::of)
        .unwrap_or_default();
    let sort = view.report_sort(formatting.formatting.report_sort());

    // Each query's Seq-aligns as the report reads them (`PrepareBlastUngappedSeqalign` where the
    // search is ungapped), their subjects and the runs of their subject ids, up to a Seq-align
    // whose Seq-ids do not read (`ids_stop`): the report ends in that query's, after the queries
    // before it.
    let mut runs = RunSubjects::default();
    let mut planned: Vec<Vec<(ReadAlign<'_>, u32)>> = Vec::with_capacity(sets.len());
    let mut stopped: Option<(usize, anyhow::Error)> = None;
    'queries: for (q, set) in sets.iter().enumerate() {
        let all = read_set(&alignments[set.clone()], ungapped && !set.is_empty());
        let mut aligns = Vec::with_capacity(all.len());
        let mut previous: Option<(&ParsedSeqId, usize)> = None;
        // The runs of subject ids among the query's Seq-aligns before this one.
        let mut runs_before = 0usize;
        for align in all {
            if matches!(align.variant(), "std" | "disc") {
                return Err(refused(
                    "Std-seg or discontinuous Seq-aligns (those of a translated search)",
                ));
            }
            let read_by_table =
                sort.sorts() || runs_before.saturating_sub(1) < formatting.counts.num_descriptions;
            let query_id = match align.seq_id(0) {
                Ok(id) => id,
                Err(stop) => {
                    stopped = Some((q, ids_stop(stop, read_by_table)?));
                    break 'queries;
                }
            };
            if scope.queries_holding(query_id).as_slice() != [q] {
                return Err(refused(
                    "a query id that names another query Bioseq of the request, or none, or \
                     several",
                ));
            }
            let id = match align.seq_id(1) {
                Ok(id) => id,
                Err(stop) => {
                    stopped = Some((q, ids_stop(stop, read_by_table)?));
                    break 'queries;
                }
            };
            let source = subjects.index_of(inputs, &scope, id)?;
            let s_idx = match previous {
                Some((previous_id, _)) if compare(previous_id, id) == Sic::Yes => {
                    runs.sources.len() - 1
                }
                Some((_, previous_source)) if previous_source == source => {
                    return Err(refused(
                        "consecutive Seq-aligns of one subject under different ids, which NCBI \
                         BLAST+ shows as two subjects",
                    ))
                }
                _ => {
                    runs.sources.push(source);
                    runs_before += 1;
                    runs.sources.len() - 1
                }
            };
            previous = Some((id, source));
            aligns.push((align, s_idx as u32));
        }
        planned.push(aligns);
    }

    // The values of each source as the writers show them, checked as the search programs check
    // the subjects that their reports show.
    let records = BioseqRecords::of(
        &subjects.sources,
        &scope,
        inputs.believe,
        inputs.database.is_some(),
        false,
        query_sets_taxon(&request.queries),
    )?;
    let db = inputs.database.map(|(_, db)| db);
    let taxdb_found = db.is_some_and(|db| inputs.taxdb.found(db));
    let mut shown: Vec<SourceShown> = Vec::with_capacity(subjects.sources.len());
    for (at, source) in subjects.sources.iter().enumerate() {
        shown.push(source_shown(
            source,
            at,
            &records,
            &scope,
            db,
            protein_subjects,
            taxdb_found,
            formatting.show_gis,
        )?);
    }
    let run_parsed: Vec<Option<Arc<ParsedSubject>>> = runs
        .sources
        .iter()
        .map(|&at| shown[at].parsed.clone())
        .collect();
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:913-917
    // Every subject of a query may print a description line (a superset of those that the
    // counts show: the refusal is the same).
    let reported = |q: usize| -> Vec<u32> {
        let mut found: Vec<u32> = planned[q].iter().map(|&(_, s_idx)| s_idx).collect();
        found.dedup();
        found
    };
    match db {
        Some(db) => refuse_described_lines_of(
            planned.len(),
            reported,
            |q| Ok(reported(q)),
            &run_parsed,
            |_| HitLoad::Loaded,
            |s_idx| shown[runs.sources[s_idx as usize]].oid.unwrap_or(db.oid(0)),
            protein_subjects,
            PROGRAM,
        )?,
        None => refuse_described_records_of(
            planned.len(),
            reported,
            |q| Ok(reported(q)),
            &run_parsed,
            protein_subjects,
            PROGRAM,
        )?,
    }

    // The texts of the loads of the database subjects at the prefetch (`tabular`'s rule).
    let load_texts = match db {
        Some(db) => prefetch_messages(
            db,
            inputs.taxdb,
            subjects.sources.iter().filter_map(|source| match source {
                SubjectSource::Database {
                    ordinal, header, ..
                } => Some((
                    *ordinal,
                    match header {
                        ReportedHeader::Parsed(subject) => Some(&**subject),
                        ReportedHeader::Stored(_) => None,
                    },
                )),
                SubjectSource::Bioseq(_) | SubjectSource::Query(_) => None,
            }),
            PROGRAM,
        )?,
        None => Vec::new(),
    };
    if !load_texts.is_empty() && (!prefetched || subjects.sources.len() != prefetched_subjects) {
        return Err(refused(
            "texts of the taxonomy database at the loads of the subjects, which the report \
             makes where the prefetch does not load them (PRE_FETCH_SEQS_LIMIT, or Seq-aligns \
             that NCBI BLAST+'s prefetch cannot read)",
        ));
    }

    // The HSPs as `CDisplaySeqalign` shows them, before any byte.
    let query_masks = query_masks(&queries, sets, pairwise.masks, &planned)?;
    let mut letters = LetterCache::default();
    let mut map = MapBuffers::default();
    let mut hits: Vec<PairwiseHit> = Vec::with_capacity(seq_aligns);
    // The subject letters that each run's rows read, and the query of the run.
    let mut spans: Vec<(usize, Option<(u32, u32)>)> = vec![(0, None); runs.sources.len()];
    for (q, aligns) in planned.iter().enumerate() {
        for (align, s_idx) in aligns {
            let source_at = runs.sources[*s_idx as usize];
            let source = &subjects.sources[source_at];
            let ds = align.dense_seg().map_err(read_stop)?;
            let range = subject_range(&ds);
            let span = &mut spans[*s_idx as usize];
            span.0 = q;
            span.1 = match (span.1, range) {
                (Some((from, to)), Some((a, b))) => Some((from.min(a), to.max(b))),
                (span, range) => span.or(range),
            };
            let row_letters = match letters.letters(range, &scope, inputs, source, q) {
                Ok(letters) => letters,
                Err(LetterEnd::Refusal(error)) => return Err(error),
                Err(LetterEnd::Exception(_)) => {
                    return Err(refused(
                        "a sequence whose letters NCBI BLAST+'s alignment display cannot read",
                    ))
                }
            };
            let context = HitContext {
                nucleotide,
                matrix,
                query_length: usize::try_from(queries[q].length).unwrap_or(0),
                subject_length: shown[source_at].length,
                title: &shown[source_at].title,
            };
            let hit = pairwise_hit(
                align,
                &ds,
                q as u32,
                *s_idx,
                &context,
                &row_letters,
                &mut map,
            )
            .map_err(read_stop)?;
            hits.push(hit);
        }
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4290-4320
    // The masks lower the query rows of outfmt 0 (`x_FillLocList`); the query-anchored reports
    // read them with the sequences (`ArchiveSequences`).
    let html = formatting.html;
    let mut html_masks = None;
    let sequences = if view == PairwiseView::Pairwise {
        html_masks = html.then(Vec::new);
        crate::algorithm::blastp::blast_engine::set_row_masks(
            &mut hits,
            Some(query_masks.as_slice()),
            [&HashMap::new(), &HashMap::new()],
            html_masks.as_mut(),
        );
        None
    } else {
        Some(ArchiveSequences::of(
            &queries,
            &subject_bioseqs,
            &subjects.sources,
            &runs,
            &shown,
            &spans,
            &mut letters,
            &scope,
            inputs,
            &query_masks,
            formatting.show_gis,
        )?)
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:90, 224
    let links = crate::blastinput::app::html_links_of(
        html,
        "blast_formatter",
        program_type.subject_nucleotide,
        inputs.database.map(|(name, _)| name),
        &[inputs.choice.number],
    )?;
    let config = PairwiseConfig {
        line_length: formatting.formatting.line_length(),
        show_frame: false,
        program: program.to_string(),
        protein_matrix: matrix,
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:423-424
        sort,
        view,
        // The OIDs of a database object's subjects by subject index (a query Bioseq's self hit
        // has its own ids, which name it); none for `-subject` Bioseqs.
        html: links.as_ref().map(|links| {
            let shown = &shown;
            crate::report::html_ids::HtmlReport::new(
                links,
                db.into_iter().flat_map(|db| {
                    runs.sources
                        .iter()
                        .map(move |&at| shown[at].oid.unwrap_or(db.oid(0)))
                }),
            )
            // A Bioseq's local id `Subject_N` (NCBI's anchor and id link) is the search's.
            .with_record_ids(
                runs.sources
                    .iter()
                    .map(|&at| {
                        shown[at]
                            .parsed
                            .as_deref()
                            .filter(|subject| subject.shown_as_record())
                            .map(|subject| subject.alignment_id().clone())
                    })
                    .collect(),
            )
        }),
        unreadable_query: None,
    };
    let subject_titles: Vec<Arc<[u8]>> = runs
        .sources
        .iter()
        .map(|&at| Arc::clone(&shown[at].title))
        .collect();
    let subject_sources: Vec<SourceTaxname> = runs
        .sources
        .iter()
        .map(|&at| shown[at].source.clone())
        .collect();
    let query_titles: Vec<Arc<[u8]>> = queries
        .iter()
        .map(|query| Arc::from(query.title.as_slice()))
        .collect();
    if let Some(sequences) = &sequences {
        refuse_anchored_aborts(
            view,
            formatting.formatting.report_sort(),
            formatting.counts.num_alignments,
            sequences,
            &hits,
            &subject_titles,
            nucleotide,
            ungapped,
        )?;
    }
    let anchored = sequences
        .as_ref()
        .map(|sequences| sequences as &dyn ReportSequences);

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2187-2250
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:126-128
    let blocks = results.karlin_blocks();
    let stats = results.search_stats(&request.service);
    let ancillary = Ancillary::of_archive(
        blocks.ungapped,
        blocks.gapped,
        stats.effective_search_space,
        psi_task,
        matches!(options.program, EProgram::PsiBlast | EProgram::DeltaBlast),
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:131-154
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:533-552
    // The database as it is now, or the `-subject` Bioseqs (`FillScanModeBlastDbInfo`: their
    // number and their total length cast to an `int`, no name).
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:363-370,402-416
    // The formatter's database lines carry the archive's filtering algorithm (`Masked using:`),
    // which `opened_database` found in every name's database.
    let filtering = options.db_filtering_algorithm_id.filter(|&id| id != -1);
    let database: Option<PairwiseDatabase> = match (db, filtering) {
        (Some(db), Some(id)) => {
            let shown = db
                .name_mask_details(id)
                .map_err(|_| refused("a filtering algorithm that the database does not hold"))?;
            Some(pairwise_database_masked(db, Some(&shown)))
        }
        (db, _) => db.map(pairwise_database),
    };
    let (num_sequences, total_letters) = {
        let total: i64 = subject_bioseqs
            .iter()
            .map(|bioseq| bioseq.length().unwrap_or(0))
            .fold(0i64, i64::wrapping_add);
        (subject_bioseqs.len() as i64, i64::from(total as i32))
    };
    let counts = &formatting.counts;

    let cell = SinkCell {
        sink: RefCell::new(sink),
        failed: RefCell::new(None),
    };
    let mut out = StreamOut(&cell);
    let mut posts = PostOut(&cell);
    let mut warnings = QueryWarnings::of_posts(&[], None, &mut posts);
    let loads = HitLoads::default();
    let row_masks = html_masks.as_deref();
    let gapped_karlin = ancillary.gapped.unwrap_or(KarlinParams {
        lambda: 0.0,
        k: 0.0,
        h: 0.0,
        alpha: 0.0,
        beta: 0.0,
    });
    let mut blastn_report = nucleotide.then(|| BlastnPairwiseReport {
        version: crate::wording::REPORT_VERSION.to_string(),
        intro: BlastnIntro {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:417-418
            megablast: options.program == EProgram::Megablast,
            indexed_megablast: false,
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:168-181
            matrix_only: values.reward == 0 && values.penalty == 0,
        },
        database_name: scan_mode_definition_bytes(b""),
        database_num_sequences: num_sequences,
        database_total_letters: total_letters,
        database: database.clone(),
        reward: values.reward,
        penalty: values.penalty,
        gap_open: values.gap_open,
        gap_extend: values.gap_extend,
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2270-2275
        // (`m_Program` is "blastn" for every task of a nucleotide search.)
        zero_gap_extension_formula: true,
        window_size: values.window_size,
        gapped: values.gapped,
        num_descriptions: counts.num_descriptions,
        num_alignments: counts.num_alignments,
        unsearched: vec![false; queries.len()],
        epilog: false,
        prolog: true,
        parsed_subjects: run_parsed.clone(),
        subject_sources: subject_sources.clone(),
        subject_oids: db.map(|db| {
            runs.sources
                .iter()
                .map(|&at| shown[at].oid.unwrap_or(db.oid(0)).get())
                .collect()
        }),
    });
    let blastn_queries: Vec<BlastnPairwiseQuery> = queries
        .iter()
        .map(|query| BlastnPairwiseQuery {
            query_name: query.title.clone(),
            query_length: usize::try_from(query.length).unwrap_or(0),
            karlin: None,
            effective_search_space: 0,
            ancillary: Some(ancillary),
        })
        .collect();
    let blastp_report = BlastpPairwiseReport {
        version: crate::wording::REPORT_VERSION.to_string(),
        database_name: scan_mode_definition_bytes(b""),
        database_num_sequences: num_sequences,
        database_total_letters: total_letters,
        matrix_name: pairwise
            .epilog_matrix
            .map(|name| String::from_utf8_lossy(name).into_owned()),
        gap_open: values.gap_open,
        gap_extend: values.gap_extend,
        word_threshold: f64::from(values.word_threshold),
        window_size: values.window_size,
        gapped_karlin,
        gumbel: None,
        gapped: values.gapped,
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:164-166
        show_linked_set_size: values.sum_statistics && ungapped,
        num_descriptions: counts.num_descriptions,
        num_alignments: counts.num_alignments,
    };
    let blastp_queries: Vec<BlastpPairwiseQuery> = queries
        .iter()
        .map(|query| BlastpPairwiseQuery {
            query_name: Vec::new(),
            query_length: usize::try_from(query.length).unwrap_or(0),
            ungapped_karlin: gapped_karlin,
            effective_search_space: 0,
            valid: true,
            batch_skipped: false,
            ancillary: Some(ancillary),
        })
        .collect();
    let round_number = match inputs.iteration {
        Ok(Some(number)) => number,
        _ => 0,
    };
    // One call of the program's writer: the prolog where `prolog`, the queries `range`, the
    // epilog where `epilog`.
    let mut write = |range: Range<usize>, prolog: bool, epilog: bool| -> io::Result<()> {
        if let Some(report) = blastn_report.as_mut() {
            report.prolog = prolog;
            report.epilog = epilog;
            return write_blastn_pairwise_report(
                &hits,
                &mut out,
                &config,
                &blastn_queries,
                range,
                &subject_titles,
                report,
                None,
                Some(&mut warnings),
                None,
                anchored,
                row_masks,
                &loads,
            );
        }
        if psi_task {
            write_psiblast_round_report(
                &hits,
                range,
                &mut out,
                &config,
                &blastp_queries,
                &query_titles,
                &subject_titles,
                &run_parsed,
                &subject_sources,
                &loads,
                &blastp_report,
                database.as_ref(),
                None,
                Some(&mut warnings),
                prolog,
                ReportRound {
                    number: round_number,
                    found_before: None,
                },
                anchored,
            )?;
            if epilog {
                write_psiblast_epilog(
                    &mut out,
                    &blastp_report,
                    database.as_ref(),
                    config.html.is_some(),
                )?;
            }
            return Ok(());
        }
        write_blastp_pairwise_window(
            &hits,
            range,
            &mut out,
            &config,
            &blastp_queries,
            &query_titles,
            &subject_titles,
            &run_parsed,
            &subject_sources,
            &blastp_report,
            None,
            database.as_ref(),
            None,
            Some(&mut warnings),
            prolog,
            epilog,
            None,
            anchored,
            row_masks,
            &loads,
        )
    };

    // `PrintProlog`; the prefetch (its variable read, its warning and the texts of its loads);
    // each query's `PrintOneResultSet` (psiblast's round number read at the first: a
    // `Web_StepNumber` that is no integer throws there, `GetPsiNumberOfIterations`); the epilog.
    write(0..0, true, false).map_err(|error| cell.error_of(error))?;
    if let Err(error) = prefetch_limit {
        return Err(error);
    }
    {
        let mut sink = cell.sink.borrow_mut();
        if prefetches && !prefetched {
            sink.post(&prefetch_warning_of(None));
        }
        sink.post(&load_texts);
    }
    if psi_task && inputs.iteration.is_err() && !sets.is_empty() {
        return Err(super::archive::invalid_archive());
    }
    let Some((q, error)) = stopped else {
        return write(0..planned.len(), false, true).map_err(|error| cell.error_of(error));
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1474-1494
    // The queries before, then the start of the report of the query whose Seq-align throws.
    write(0..q, false, false).map_err(|error| cell.error_of(error))?;
    crate::report::pairwise::write_query_preamble(
        &mut out,
        psi_task.then_some(round_number),
        &queries[q].title,
        usize::try_from(queries[q].length).unwrap_or(0),
        config.html.is_some(),
    )
    .and_then(|()| out.flush())
    .map_err(|error| cell.error_of(error))?;
    Err(error)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1082-1128
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1541
/// The end of the report at a Seq-align whose Seq-ids do not read (`stop`): one of NCBI's
/// exceptions (`archive::stop_error`) where NCBI reads them after the query's `Query=` block
/// (`read_by_table`): the sort of `-sorthits` or `-sorthsps` above 0 reads every Seq-align's, the
/// description table those of its first `-num_descriptions` subjects and the first Seq-align of
/// the next. A Seq-align past them is read by the alignments' display, which is refused, as the
/// stops NCBI has no result for.
fn ids_stop(stop: Stop, read_by_table: bool) -> Result<anyhow::Error> {
    match stop {
        Stop::Exception(_) | Stop::InvalidSelection { .. } | Stop::Unassigned { .. }
            if read_by_table =>
        {
            Ok(super::archive::stop_error(stop))
        }
        Stop::Exception(_) | Stop::InvalidSelection { .. } | Stop::Unassigned { .. } => {
            Err(refused(
                "a Seq-align whose Seq-ids do not read past the Seq-aligns of the description \
                 table, which NCBI BLAST+'s alignment display reads",
            ))
        }
        other => Err(super::archive::stop_error(other)),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4545-4572
/// The sequences of an object's query-anchored reports (outfmt 1 to 4,
/// `query_anchored::ReportSequences`), read as outfmt 0's rows read them:
///
/// - Letters: each query Bioseq's whole, each run subject's from the first to the last letter
///   that its rows read (`LetterCache`), at the sequences' positions.
/// - Labels: `x_FillSeqid` of the scope's Bioseq: a query or `-subject` Bioseq by its stored
///   Seq-ids (`query_anchored::fill_seqid`), a database sequence by its parsed Seq-ids, else its
///   OID (`query_anchored::database_label`).
/// - Masks: the query's masks that its rows show (`query_masks`); the archive keeps no subject
///   masks.
struct ArchiveSequences<'a> {
    queries: Vec<Vec<u8>>,
    /// By run subject index: where its letters start, and its letters.
    subjects: Vec<(u32, Vec<u8>)>,
    query_labels: Vec<Vec<u8>>,
    subject_labels: Vec<Vec<u8>>,
    /// By run subject index: the query whose Bioseq it is (a self-hit, `SubjectSource::Query`).
    self_hits: Vec<Option<usize>>,
    query_masks: &'a [Vec<(i32, i32)>],
}

impl<'a> ArchiveSequences<'a> {
    #[allow(clippy::too_many_arguments)]
    fn of(
        queries: &[QueryEntry<'_>],
        subject_bioseqs: &[crate::report::archive::read::request::Bioseq<'_>],
        sources: &[SubjectSource],
        runs: &RunSubjects,
        shown: &[SourceShown],
        spans: &[(usize, Option<(u32, u32)>)],
        letters: &mut LetterCache,
        scope: &ScopeBioseqs<'_, '_>,
        inputs: &ObjectInputs<'_>,
        query_masks: &'a [Vec<(i32, i32)>],
        show_gis: bool,
    ) -> Result<Self> {
        let label = |ids: &mut dyn Iterator<Item = &ParsedSeqId>| {
            let ids: Vec<ParsedSeqId> = ids.cloned().collect();
            let mut out = Vec::new();
            fill_seqid(&ids, show_gis, &mut out);
            out
        };
        let unread = |end: LetterEnd| match end {
            LetterEnd::Refusal(error) => error,
            LetterEnd::Exception(_) => {
                refused("a sequence whose letters NCBI BLAST+'s alignment display cannot read")
            }
        };
        let mut query_letters = Vec::with_capacity(queries.len());
        for query in queries {
            query_letters
                .push(super::tabular::read_bioseq(query.bioseq, "a query Bioseq").map_err(unread)?);
        }
        let mut subjects = Vec::with_capacity(runs.sources.len());
        let mut subject_labels = Vec::with_capacity(runs.sources.len());
        let mut self_hits = Vec::with_capacity(runs.sources.len());
        for (s_idx, &at) in runs.sources.iter().enumerate() {
            let (q, span) = spans[s_idx];
            let read = letters
                .letters(span, scope, inputs, &sources[at], q)
                .map_err(unread)?;
            subjects.push((read.subject_from, read.subject.to_vec()));
            let mut out = Vec::new();
            match sources[at] {
                SubjectSource::Database { .. } => crate::report::query_anchored::database_label(
                    shown[at].parsed.as_deref(),
                    shown[at].oid.expect("a database subject's OID"),
                    &mut out,
                ),
                SubjectSource::Bioseq(bioseq) => {
                    out = label(&mut subject_bioseqs[bioseq].ids());
                }
                SubjectSource::Query(query) => {
                    out = label(&mut queries[query].bioseq.ids());
                }
            }
            subject_labels.push(out);
            self_hits.push(match sources[at] {
                SubjectSource::Query(query) => Some(query),
                _ => None,
            });
        }
        Ok(ArchiveSequences {
            queries: query_letters,
            subjects,
            query_labels: queries
                .iter()
                .map(|query| label(&mut query.bioseq.ids()))
                .collect(),
            subject_labels,
            self_hits,
            query_masks,
        })
    }
}

/// Appends the letters `from..=to` of `letters`, which start at position `start`, as far as they
/// go.
fn push_letters(letters: &[u8], start: u32, from: u32, to: u32, out: &mut Vec<u8>) {
    let (Some(from), Some(to)) = (from.checked_sub(start), to.checked_sub(start)) else {
        return;
    };
    let end = (to as usize + 1).min(letters.len());
    if let Some(range) = letters.get(from as usize..end) {
        out.extend_from_slice(range);
    }
}

impl ReportSequences for ArchiveSequences<'_> {
    fn subject_is_query(&self, q_idx: usize, s_idx: u32) -> bool {
        self.self_hits.get(s_idx as usize).copied().flatten() == Some(q_idx)
    }

    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>) {
        push_letters(&self.queries[q_idx], 0, from, to, out);
    }

    fn subject_letters(&self, s_idx: u32, from: u32, to: u32, out: &mut Vec<u8>) {
        let (start, letters) = &self.subjects[s_idx as usize];
        push_letters(letters, *start, from, to, out);
    }

    fn query_label(&self, q_idx: usize, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.query_labels[q_idx]);
    }

    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.subject_labels[s_idx as usize]);
    }

    fn query_masks(&self, q_idx: usize, out: &mut Vec<Mask>) {
        if let Some(masks) = self.query_masks.get(q_idx) {
            out.extend(
                masks
                    .iter()
                    .map(|&(from, to)| (from.max(0) as u32, to.max(0) as u32, 0)),
            );
        }
    }

    fn subject_masks(&self, _q_idx: usize, _s_idx: u32, _out: &mut Vec<Mask>) {}
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2003-2040
/// The refusal of a query-anchored report (`view`) on which NCBI BLAST+ aborts or crashes, before
/// the object's first byte, as the search programs refuse it (`pairwise::section_aborts`,
/// `pairwise::blastp_section_frees`): a section with a self-hit, whose merge leaves out the
/// subject's rows, and an ungapped protein section whose HSP of length 0 the merge frees and
/// reads.
#[allow(clippy::too_many_arguments)]
fn refuse_anchored_aborts(
    view: PairwiseView,
    sort: crate::report::hit_sort::ReportSort,
    num_alignments: usize,
    sequences: &ArchiveSequences<'_>,
    hits: &[PairwiseHit],
    subject_titles: &[Arc<[u8]>],
    nucleotide: bool,
    ungapped: bool,
) -> Result<()> {
    let io_error = |error: io::Error| anyhow::Error::from(error);
    for q in 0..sequences.queries.len() {
        let is_query = |s_idx: u32| sequences.subject_is_query(q, s_idx);
        let aborts = if !sequences.self_hits.contains(&Some(q)) {
            false
        } else if nucleotide {
            section_aborts(
                view,
                sort,
                num_alignments,
                q,
                hits,
                subject_titles,
                is_query,
                |group| blastn_section_hsps(num_alignments, group, &[], is_query),
            )
            .map_err(io_error)?
        } else {
            blastp_section_aborts(
                view,
                sort,
                num_alignments,
                q,
                hits,
                subject_titles,
                is_query,
            )
            .map_err(io_error)?
        };
        if aborts {
            return Err(refused(
                "a query Bioseq among its own subjects in the query-anchored output formats, on \
                 whose alignments NCBI BLAST+ aborts",
            ));
        }
        if ungapped
            && !nucleotide
            && blastp_section_frees(
                view,
                sort,
                num_alignments,
                q,
                hits,
                subject_titles,
                |_, _| false,
            )
            .map_err(io_error)?
        {
            return Err(refused(
                "an ungapped HSP of length 0 in the query-anchored output formats, whose merge \
                 NCBI BLAST+ frees and reads",
            ));
        }
    }
    Ok(())
}

/// The masks of each query that its rows show (`TMaskedQueryRegions` of its result,
/// `Results::masks`): those whose id is the query's Seq-align id (`x_FillLocList` compares each
/// mask's id with each row's), as 0-based inclusive ranges. A mask that names a subject of the
/// query's Seq-aligns would lower that subject's rows: refused.
fn query_masks(
    queries: &[QueryEntry<'_>],
    sets: &[Range<usize>],
    masks: &[Vec<MaskRange<'_>>],
    planned: &[Vec<(ReadAlign<'_>, u32)>],
) -> Result<Vec<Vec<(i32, i32)>>> {
    let mut out = vec![Vec::new(); queries.len()];
    for (q, aligns) in planned.iter().enumerate() {
        let Some(list) = masks.get(q).filter(|list| !list.is_empty()) else {
            continue;
        };
        if sets[q].is_empty() || aligns.is_empty() {
            continue;
        }
        let query_id = aligns[0].0.seq_id(0).map_err(super::archive::stop_error)?;
        for (align, _) in aligns {
            let subject_id = align.seq_id(1).map_err(super::archive::stop_error)?;
            if list
                .iter()
                .any(|mask| compare(mask.id, subject_id) == Sic::Yes)
            {
                return Err(refused(
                    "a query mask that names a subject of the query's alignments, whose rows \
                     NCBI BLAST+ lowers too",
                ));
            }
        }
        for mask in list {
            if compare(mask.id, query_id) != Sic::Yes {
                continue;
            }
            let (Ok(from), Ok(to)) = (i32::try_from(mask.from), i32::try_from(mask.to)) else {
                return Err(refused("a query mask past the int range"));
            };
            out[q].push((from, to));
        }
    }
    Ok(out)
}

/// The values of the subject `source` (at `at` of the object's sources) as the writers show it:
/// a database sequence's header as the search programs read a reported sequence's
/// (`ReportedHeader`, with their refusals of what outfmt 0 shows), a Bioseq's record
/// (`BioseqRecords`).
fn source_shown(
    source: &SubjectSource,
    at: usize,
    records: &BioseqRecords,
    scope: &ScopeBioseqs<'_, '_>,
    db: Option<&crate::seqdb::SeqDb>,
    protein: bool,
    taxdb_found: bool,
    show_gis: bool,
) -> Result<SourceShown> {
    let shown_parsed = |subject: &Arc<ParsedSubject>| {
        Arc::new(ParsedSubject::clone(subject).showing_gis(show_gis))
    };
    match source {
        SubjectSource::Database {
            ordinal, header, ..
        } => {
            let db = db.expect("a database subject");
            let ordinal = *ordinal;
            let oid = db.oid(ordinal);
            let length = if db.seqtype() == SeqType::Protein {
                db.protein_length(ordinal)
                    .map(|length| length as usize)
                    .map_err(|error| {
                        crate::algorithm::common::database::database_error(error, PROGRAM)
                    })?
            } else {
                db.nucleotide_length(ordinal).map_err(|error| {
                    crate::algorithm::common::database::database_error(error, PROGRAM)
                })?
            };
            match header {
                ReportedHeader::Stored(stored) => {
                    check_reported_title(
                        db,
                        ordinal,
                        stored.replaced_bytes,
                        TitleReads::Reports,
                        PROGRAM,
                    )?;
                    stored.check_shown_source(|| oid, PROGRAM)?;
                    if protein {
                        check_shown_protein_title(&stored.title, || oid, PROGRAM)?;
                    }
                    Ok(SourceShown {
                        title: Arc::from(stored.title.as_slice()),
                        parsed: None,
                        source: stored.source.clone(),
                        length,
                        oid: Some(oid),
                    })
                }
                ReportedHeader::Parsed(subject) => {
                    check_reported_title(
                        db,
                        ordinal,
                        subject.replaced_title_bytes(),
                        TitleReads::Reports,
                        PROGRAM,
                    )?;
                    check_pairwise_parsed(
                        subject,
                        HitLoad::Loaded,
                        protein,
                        taxdb_found,
                        || oid,
                        PROGRAM,
                    )?;
                    Ok(SourceShown {
                        title: Arc::from(subject.first_title()),
                        parsed: Some(shown_parsed(subject)),
                        source: SourceTaxname::default(),
                        length,
                        oid: Some(oid),
                    })
                }
            }
        }
        SubjectSource::Bioseq(_) | SubjectSource::Query(_) => {
            let title = records.records[at].title.as_slice();
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:891-905
            // A Bioseq of one local id (whose tabular rows show it as a record,
            // `BioseqRecords::ids`) shows the id before its title in the description line and
            // the heading unless it is `lcl|Subject_N`, as a `-subject` record of
            // `-parse_deflines` in the search programs (`ParsedSubject::of_record`).
            let subject = match &records.parsed[at] {
                Some(subject) => shown_parsed(subject),
                None => {
                    let bioseq = match *source {
                        SubjectSource::Bioseq(index) => scope.subjects[index],
                        SubjectSource::Query(index) => scope.queries[index].bioseq,
                        SubjectSource::Database { .. } => unreachable!("a Bioseq subject"),
                    };
                    let ids = crate::blastinput::fasta_reader::believed::BelievedIds::new(
                        bioseq.ids().cloned().collect(),
                        crate::report::archive::bioseq::SeqLocStrand::Unknown,
                    );
                    shown_parsed(&Arc::new(ParsedSubject::of_record(&ids, title)))
                }
            };
            match subject.shown_as_record() {
                false => crate::report::db_subject::check_pairwise_records(
                    std::slice::from_ref(&Some(Arc::clone(&subject))),
                    protein,
                    PROGRAM,
                )?,
                true => {
                    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:4092-4095
                    let reads_past_end = protein
                        && [false, true].into_iter().any(|leave_prefix| {
                            crate::report::defline::generate_defline(title, true, leave_prefix)
                                .is_ok_and(|defline| defline.reads_past_end)
                        });
                    if reads_past_end {
                        return Err(refused(
                            "a subject Bioseq whose title of punctuation NCBI BLAST+'s \
                             x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 \
                             crashes in outfmt 0)",
                        ));
                    }
                }
            }
            Ok(SourceShown {
                title: Arc::from(title),
                parsed: Some(subject),
                source: SourceTaxname::default(),
                length: records.lengths[at],
                oid: None,
            })
        }
    }
}

/// The subject letters that the rows of `ds` read: from the first to the last letter of row 1's
/// segments (plus strand, inclusive).
fn subject_range(ds: &PairDenseg) -> Option<(u32, u32)> {
    let mut range: Option<(u32, u32)> = None;
    for (seg, &len) in ds.lens.iter().enumerate() {
        let Ok(start) = u32::try_from(ds.starts[seg * 2 + 1]) else {
            continue;
        };
        let last = start.saturating_add(len.saturating_sub(1));
        range = Some(match range {
            Some((from, to)) => (from.min(start), to.max(last)),
            None => (start, last),
        });
    }
    range
}

/// What the HSP of a Seq-align reads besides it.
struct HitContext<'c> {
    /// A blastn report (rows compared without the matrix, the Strand line).
    nucleotide: bool,
    /// The display's matrix (`CDisplaySeqalign`'s, BLOSUM62 for an unknown or NULL name).
    matrix: ScoringMatrix,
    /// The query Bioseq's length.
    query_length: usize,
    /// The subject's length.
    subject_length: usize,
    /// The subject's title (`PairwiseHit::subject_title`).
    title: &'c [u8],
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:811-848
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1830-1862
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2120-2155
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3593-3604
/// The HSP of the Seq-align `align` (its Dense-seg `ds`) of the query `q_idx` and the report's
/// subject `s_idx` as `CDisplaySeqalign` shows it (`PairwiseHit`): its scores (`GetAlnScores`:
/// `sum_n`, `comp_adjustment_method`, `num_ident`), the rows of its Dense-seg as `CAlnVec` reads
/// them (reversed where the first strand is minus: the query on its plus strand), the counts of
/// the rows (identities on the letters as they are, positives with the display's matrix), and the
/// edit script of its segments, which the writers draw the rows' coordinates from: a Dense-seg
/// whose segments are not the script's (`seqalign::nucleotide_row_segments`,
/// `protein_row_segments`, which BLAST's Seq-aligns are) is refused.
fn pairwise_hit(
    align: &ReadAlign<'_>,
    ds: &PairDenseg,
    q_idx: u32,
    s_idx: u32,
    context: &HitContext<'_>,
    letters: &RowLetters<'_>,
    map: &mut MapBuffers,
) -> std::result::Result<PairwiseHit, ReadStop> {
    let scores = align.scores()?;
    if scores.use_this_gi {
        return Err(ReadStop::Unsupported(
            "a Seq-align with a use_this_gi score, whose GIs NCBI BLAST+ chooses the definition \
             lines by",
        ));
    }
    let minus = ds.strands.contains(&Strand::Minus);
    if ds.strands == [Strand::Minus, Strand::Minus] || (minus && !context.nucleotide) {
        return Err(ReadStop::Unsupported(
            "a Seq-align whose strands BLAST does not write",
        ));
    }
    let raw_score = i32::try_from(scores.score)
        .map_err(|_| ReadStop::Unsupported("a score past the int range"))?;
    let sum_n = match scores.sum_n {
        n if n > 0 => Some(
            i32::try_from(n).map_err(|_| ReadStop::Unsupported("a sum_n past the int range"))?,
        ),
        _ => None,
    };
    let denseg = Denseg::pair(
        [SeqKey::Query, SeqKey::Subject(0)],
        &ds.starts,
        &ds.lens,
        ds.strands,
        ds.first_minus,
    );
    let aln = AlnMap::new_in(&denseg, map).map_err(ReadStop::Unsupported)?;
    let (length, gaps, gap_opens) = aln.align_lengths();
    let row = |row: usize| -> std::result::Result<Vec<u8>, ReadStop> {
        aln.whole_aln_seq_string(row, 0, letters)
            .map(|string| string.buffer)
            .map_err(|_| ReadStop::Unsupported("a translated alignment row"))
    };
    let (query_row, subject_row) = (row(0)?, row(1)?);
    let q_start = aln.seq_start(0) as usize + 1;
    let q_end = aln.seq_stop(0) as usize + 1;
    let (s_first, s_last) = (aln.seq_start(1) as usize + 1, aln.seq_stop(1) as usize + 1);
    let (s_start, s_end) = if minus {
        (s_last, s_first)
    } else {
        (s_first, s_last)
    };
    let length = length as usize;
    let gaps = gaps as usize;
    // The edit script of the Seq-align's segments, in their order (BLAST's: the query context's
    // for a minus-strand HSP).
    let mut ops: Vec<GapEditOp> = Vec::with_capacity(ds.lens.len());
    for (seg, &len) in ds.lens.iter().enumerate() {
        let op = match (ds.starts[seg * 2] >= 0, ds.starts[seg * 2 + 1] >= 0) {
            (true, true) => GapEditOp::Sub(len),
            (false, true) => GapEditOp::Del(len),
            (true, false) => GapEditOp::Ins(len),
            (false, false) => {
                return Err(ReadStop::Unsupported(
                    "a Dense-seg segment without letters in either row",
                ))
            }
        };
        ops.push(op);
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2127-2140
    // The identities of the rows' letters as the Bioseqs give them.
    let num_ident = query_row
        .iter()
        .zip(&subject_row)
        .filter(|(query, subject)| query == subject)
        .count();
    let query_case = if query_row.iter().any(u8::is_ascii_lowercase) {
        crate::report::row_letters::QueryCase::Bioseq
    } else {
        crate::report::row_letters::QueryCase::Masked
    };
    let positives = (!context.nucleotide).then(|| {
        crate::report::row_letters::display_positives(
            &query_row,
            &subject_row,
            context.matrix,
            query_case,
        )
    });
    let hit = Hit {
        identity: if length > 0 {
            num_ident as f64 * 100.0 / length as f64
        } else {
            0.0
        },
        length,
        mismatch: length as i64 - (num_ident + gaps) as i64,
        gapopen: gap_opens as usize,
        q_start,
        q_end,
        s_start,
        s_end,
        e_value: scores.evalue,
        bit_score: scores.bit_score,
        num_ident,
        query_frame: match (context.nucleotide, minus) {
            (true, true) => -1,
            (true, false) => 1,
            (false, _) => 0,
        },
        query_length: context.query_length,
        q_idx,
        s_idx,
        raw_score,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(ops),
        num_positives: positives.unwrap_or(num_ident),
    };
    // The segments that the writers make of the HSP are the Seq-align's.
    let mut segments = Segments::default();
    if context.nucleotide {
        crate::report::serial::seqalign::nucleotide_row_segments(&hit, &mut segments);
    } else {
        crate::report::serial::seqalign::protein_row_segments(&hit, &mut segments);
    }
    let same = segments.starts() == ds.starts.as_slice()
        && segments
            .lens()
            .iter()
            .map(|&len| len as u32)
            .eq(ds.lens.iter().copied());
    if !same {
        return Err(ReadStop::Unsupported(
            "a Dense-seg whose segments are not those of an edit script of BLAST's",
        ));
    }
    let [query_stored_gaps, subject_stored_gaps] = if context.nucleotide {
        [Vec::new(), Vec::new()]
    } else {
        [
            stored_gap_columns(ds, &query_row, 0),
            stored_gap_columns(ds, &subject_row, 1),
        ]
    };
    Ok(PairwiseHit {
        hit,
        query_seq: Some(String::from_utf8_lossy(&query_row).into_owned()),
        subject_seq: Some(String::from_utf8_lossy(&subject_row).into_owned()),
        query_frame: (!context.nucleotide).then_some(1),
        subject_frame: (!context.nucleotide).then_some(1),
        positives,
        gaps: Some(gaps),
        subject_length: Some(context.subject_length),
        subject_title: crate::blastinput::fasta_reader::hit_subject_title(context.title),
        comp_adjust_method: u8::try_from(scores.comp_adjustment_method).ok(),
        sum_n,
        subject_stored_gaps,
        query_stored_gaps,
        query_case,
        seqalign_num_ident: scores.num_ident,
        edge: PairwiseEdge::Normal,
    })
}

/// The columns (0-based) of the alignment of `ds` where the row `letters` of `row` shows a
/// residue that is the gap character (`-`: a PSSM query's gap residue, a database protein's
/// stored gap): a `-` in a segment where the row has letters.
fn stored_gap_columns(ds: &PairDenseg, letters: &[u8], row: usize) -> Vec<u32> {
    if !letters.contains(&b'-') {
        return Vec::new();
    }
    // The display's segments run in the order of the reversed Dense-seg where its first strand
    // is minus (a protein row has none).
    let mut columns = Vec::new();
    let mut column = 0usize;
    for (seg, &len) in ds.lens.iter().enumerate() {
        let len = len as usize;
        if ds.starts[seg * 2 + row] >= 0 {
            for at in column..column + len {
                if letters.get(at) == Some(&b'-') {
                    columns.push(at as u32);
                }
            }
        }
        column += len;
    }
    columns
}

/// The report's sink shared by the writer's stream (`StreamOut`) and its posts (`PostOut`), with
/// the error of the stream where a write failed.
struct SinkCell<'s> {
    sink: RefCell<&'s mut dyn ReportSink>,
    failed: RefCell<Option<anyhow::Error>>,
}

/// The error of a write that the sink's stream refused.
const STREAM_FAILED: &str = "blast_formatter's output stream failed";

impl SinkCell<'_> {
    fn keep(&self, result: Result<()>) -> io::Result<()> {
        result.map_err(|error| {
            *self.failed.borrow_mut() = Some(error);
            io::Error::other(STREAM_FAILED)
        })
    }

    /// The run's error of a writer's `error`: the stream's where it failed.
    fn error_of(&self, error: io::Error) -> anyhow::Error {
        match self.failed.borrow_mut().take() {
            Some(failed) => failed,
            None => crate::cli::NativeError::of_writer(error),
        }
    }
}

/// The `-out` stream as the writers' `ReportWrite`: the bytes inserted, NCBI's flushes, the
/// stream's `badbit`, each failure the stream's error (`BLAST failed to write output`).
struct StreamOut<'c, 's>(&'c SinkCell<'s>);

impl Write for StreamOut<'_, '_> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let result = self.0.sink.borrow_mut().write(buf);
        self.0.keep(result)?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl ReportWrite for StreamOut<'_, '_> {
    fn ncbi_flush(&mut self) -> io::Result<()> {
        let result = self.0.sink.borrow_mut().flush();
        self.0.keep(result)
    }

    fn badbit(&mut self) -> io::Result<()> {
        let error = self.0.sink.borrow_mut().badbit();
        self.0.keep(Err(error))
    }
}

/// The report's posts (`QueryWarnings`' sink): each text posted after the stream's bytes before
/// it.
struct PostOut<'c, 's>(&'c SinkCell<'s>);

impl Write for PostOut<'_, '_> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.sink.borrow_mut().post(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
#[path = "pairwise_tests.rs"]
mod tests;
