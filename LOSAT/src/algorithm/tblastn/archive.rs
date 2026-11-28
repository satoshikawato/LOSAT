//! TBLASTN's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:294-297
//! One archive per query batch (the reader's batch, `BATCH_SIZE` letters, 20000 by default),
//! written after the batch's search. LOSAT's TBLASTN searches the batches first and writes every
//! format at the end, from the hit lists it holds (`stage_e_report::render`), so the archives are
//! written there, one per retained batch with that batch's searched flag; `report::archive` writes
//! the document and this module gathers what TBLASTN holds for its other formats: the Seq-aligns of
//! outfmt 8, the SEG and lower-case masks of the queries, and the Karlin-Altschul blocks and search
//! space of the batch's first query. A query is a protein with one context, so the first valid
//! context of a batch is its first query (a first query without one gives both blocks the `-1`
//! triple), and the batch's `EffectiveSearchSpace` option is that query's search space.

use std::ops::Range;

use anyhow::Result;

use super::seqalign::query_records_of;
use super::stage_d_results::KappaResultHitList;
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::{ReportRanges, ShownSubjects};
use super::x_subject_store::SubjectLetters;
use super::TblastnArgs;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::core::blast_util::ncbi_int4_from_double;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{OptionInputs, RequestKind};
use crate::report::archive::read::import;
use crate::report::archive::read::options::{EProgram, Setter};
use crate::report::archive::requests::{LcaseFrames, RequestPolicy, RunRequests};
use crate::report::archive::{
    ka_block, write_archive, AlignSource, ArchiveHandle, BatchArchive, KaBlock,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::serial::seqalign::{write_seq_aligns, Segments};
use crate::stats::tables::KarlinParams;
use crate::utils::seg::SegParams;

/// What a TBLASTN run's arguments give its search strategy's requests
/// (`strategy::StrategyExport::for_args`): the archive's recorded options, the counts
/// of the `HitlistSize` and the database of `-db`.
impl crate::report::archive::strategy::StrategyArgs for TblastnArgs {
    fn strategy_inputs(&self) -> crate::report::archive::strategy::StrategyInputs {
        crate::report::archive::strategy::StrategyInputs {
            inputs: option_inputs(self),
            counts: crate::report::archive::strategy::HitlistArgs {
                max_target_seqs: self.max_target_seqs,
                num_descriptions: self.num_descriptions,
                num_alignments: self.num_alignments,
            },
            database: self.database_name(),
        }
    }
}

/// The recorded options of a TBLASTN run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2). A run of `-in_pssm` is PSI-TBLASTN's
/// (`RequestKind::PsiTblastn`: psiblast's handle, whose rules `algorithm_options` applies to the
/// same arguments).
pub(crate) fn option_inputs(args: &TblastnArgs) -> OptionInputs {
    let kind = if args.in_pssm.is_some() {
        RequestKind::PsiTblastn
    } else {
        RequestKind::from(Program::Tblastn)
    };
    let mut inputs = OptionInputs::new(kind, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    inputs.evalue = args.evalue;
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.matrix = args
        .matrix
        .as_deref()
        .map(|v| crate::utils::os_bytes::os_str_bytes(v).into_owned());
    inputs.xdrop_gap = args.xdrop_gap;
    inputs.xdrop_gap_final = args.xdrop_gap_final;
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps.map(|count| count as i32);
    inputs.culling_limit = args.culling_limit;
    inputs.best_hit_overhang = args.best_hit_overhang;
    inputs.best_hit_score_edge = args.best_hit_score_edge;
    inputs.subject_besthit = args.subject_besthit;
    inputs.sum_stats = args.sum_stats;
    inputs.use_sw_tback = args.use_sw_tback;
    // `-threshold` is read as a double and recorded as an `int` (`static_cast<int>`), converted
    // as the x86-64 build does (`INT_MIN` past `Int4`, as the lookup table takes it:
    // `ncbi_int4_from_double`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:884-892
    inputs.threshold = args.threshold.map(ncbi_int4_from_double);
    inputs.window_size = args.window_size;
    inputs.ungapped = args.ungapped;
    inputs.soft_masking = args.soft_masking;
    inputs.seg = args.seg.as_deref().and_then(|v| {
        crate::report::archive::options::seg_filter(crate::utils::os_bytes::os_str_bytes(v))
    });
    inputs.max_intron_length = args.max_intron_length.unwrap_or(0);
    inputs.db_gencode = i32::from(args.db_gencode);
    inputs.imported = args.import.as_ref().map(|run| run.recorded.clone());
    inputs.application = args.import.as_ref().map(|run| run.app.ncbi_name());
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(args.formatting.recorded_hitlist_start(), |size| size as i32);
    if let Ok((mode, _)) = crate::blastinput::app::parse_comp_based_stats(
        crate::utils::os_bytes::os_str_bytes(&args.comp_based_stats),
        false,
        false,
    ) {
        inputs.comp_based_stats = mode as i32;
    }
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with TBLASTN's recorded options).
pub(crate) fn handle_for(
    args: &TblastnArgs,
    archive_requested: bool,
    options: &crate::report::archive::OptionsHandle,
) -> Option<ArchiveHandle> {
    if !archive_requested {
        return None;
    }
    // NCBI keeps the name without the white space at its ends (`NStr::TruncateSpaces`,
    // blast_args.cpp:2472-2476), and the request's `database` is that name.
    Some(ArchiveHandle::for_program(
        &option_inputs(args),
        args.database_name().as_deref(),
        options,
    ))
}

/// The requests of a run's batches (`RunRequests`): protein queries, nucleotide subjects (the
/// `-subject` records `subjects`, `None` for `-db`), and with `-lcase_masking` (`mask_lowercase`)
/// one `notset` mask of each query's lower-case runs.
pub(super) fn requests(subjects: Option<&[FastaRecord]>, mask_lowercase: bool) -> RunRequests<'_> {
    RunRequests::new(
        RequestPolicy {
            query_protein: true,
            subject_protein: false,
            lcase: mask_lowercase.then_some(LcaseFrames::NotSet),
        },
        subjects,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1718-1725
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
/// The `EffectiveSearchSpace` option of a batch (`BatchArchive::searched_space`,
/// `StrategyExport::note_searched`): the search space of its first query's context, `first`
/// (`-searchsp` in an invalid context too), when NCBI searched the batch (`searched`: a query of
/// it has a valid context), and that context's query length (`archive::RequestLength::space`):
/// the letters of the first query, `record`, or 0 where its context is not valid (`valid`).
pub(super) fn batch_searched_space(
    searched: bool,
    first: &super::stage_d_stats::LocalContextLength,
    record: &FastaRecord,
    valid: bool,
) -> Option<crate::report::archive::SearchedContext0> {
    searched.then(|| crate::report::archive::SearchedContext0 {
        space: first.eff_searchsp,
        query_length: if valid { record.seq().len() as i64 } else { 0 },
    })
}

pub(crate) use crate::report::archive::error_lines;

/// What the archives of a run need besides the report's own data.
pub(super) struct ArchiveRender<'a> {
    /// The archive (outfmt 11).
    pub handle: &'a ArchiveHandle,
    /// The searched queries of each batch NCBI read and searched, in order (the batches before
    /// a stop; `read_query_batches`).
    pub spans: &'a [Range<usize>],
}

/// The values of a batch's first query that `search-stats` and `ka-blocks` take (the
/// `EffectiveSearchSpace` option: `batch_searched_space`).
struct FirstContext {
    /// The ungapped and the gapped block; `None` when the first query has no valid context
    /// (both blocks are then missing). An ungapped search has no gapped block (`None`).
    karlin: Option<(KarlinParams, Option<KarlinParams>)>,
    /// The effective search space and length adjustment of that context (0 without one).
    search_space: i64,
    length_adjustment: i64,
}

/// The masked regions of each query as NCBI's `TMaskedQueryRegions` holds them (`Map`,
/// api/blast_aux.cpp:825-842): the SEG intervals and, with `-lcase_masking`, the lower-case runs of
/// the searched letters, merged where they overlap, a mask that covers the whole query dropped,
/// and the rest moved into the query's record by the start of its interval; inclusive 0-based
/// positions.
///
/// The masks are computed here, per batch, at report time: SEG runs again on every query of the
/// batch (and the lower-case runs are read again), instead of the search keeping each query's
/// mask until the last archive is written. That is the price of not holding the masks of all
/// batches (memory); it is CPU only, for outfmt 11 only, and below
/// NCBI's cost (200,000 peptides: 24 s user, NCBI 87 s).
fn query_record_masks(
    queries: &[FastaRecord],
    first: usize,
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    ranges: &ReportRanges<'_>,
) -> Vec<Vec<(i32, i32)>> {
    queries
        .iter()
        .enumerate()
        .map(|(offset, query)| {
            let q_idx = first + offset;
            let letters = query.seq();
            let intervals =
                super::search_seed::tblastn_report_mask_intervals(letters, seg, mask_lowercase)
                    .into_iter()
                    .map(|(start, end)| (start as i32, end as i32 - 1))
                    .collect();
            let offset = ranges.queries.offset(q_idx) as i32;
            crate::algorithm::common::masks::combine_masks(intervals)
                .into_iter()
                .filter(|&(from, to)| from <= to && !(from == 0 && to == letters.len() as i32 - 1))
                .map(|(from, to)| (from + offset, to + offset))
                .collect()
        })
        .collect()
}

/// The Seq-aligns of one batch: the elements of outfmt 8, one query at a time.
struct TblastnAligns<'a> {
    run: RunRecord<'a>,
    queries: &'a [FastaRecord],
    first: usize,
    hitlists: &'a [KappaResultHitList],
    subject_records: &'a [FastaRecord],
    subject_letters: SubjectLetters<'a>,
    shown: &'a ShownSubjects<'a>,
    ranges: &'a ReportRanges<'a>,
}

impl<'a> AlignSource for TblastnAligns<'a> {
    fn write<E: crate::report::serial::encoder::Encoder>(
        &mut self,
        enc: &mut E,
        flush: &mut dyn FnMut(&mut E),
    ) -> Result<()> {
        let mut segments = Segments::default();
        let mut subjects: Vec<SubjectRecord<'a>> = Vec::new();
        let mut hsps: Vec<HspRecord<'a>> = Vec::new();
        let hitlists: &'a [KappaResultHitList] = self.hitlists;
        for (offset, hitlist) in hitlists.iter().enumerate() {
            let q_idx = self.first + offset;
            query_records_of(
                hitlist,
                self.subject_records,
                self.subject_letters,
                self.shown,
                self.ranges,
                &mut subjects,
                &mut hsps,
            )?;
            if hsps.is_empty() {
                continue;
            }
            let query = &self.queries[q_idx];
            let query = QueryRecord {
                id: query.seq_id(),
                searched_length: i32::try_from(query.seq().len())?,
                offset: i32::try_from(self.ranges.queries.offset(q_idx))?,
            };
            write_seq_aligns(
                enc,
                &self.run,
                &query,
                &subjects,
                &hsps,
                &mut segments,
                |enc| flush(enc),
            );
        }
        Ok(())
    }
}

/// The report data of a run, as `stage_e_report::render` holds it.
pub(super) struct ArchiveData<'a> {
    pub queries: &'a [FastaRecord],
    pub subject_records: &'a [FastaRecord],
    pub subject_letters: SubjectLetters<'a>,
    pub shown: &'a ShownSubjects<'a>,
    pub hitlists: &'a [KappaResultHitList],
    pub parameters: &'a LocalSubjectParameters,
    pub ungapped_karlin: &'a [KarlinParams],
    pub query_validity: &'a [bool],
    pub query_batch_skipped: &'a [bool],
    /// The gapped Karlin-Altschul values; `None` for an ungapped search, which makes no gapped
    /// block (NCBI's `kbp_gap` is not set up).
    pub gapped: Option<KarlinParams>,
    pub seg: Option<&'a SegParams>,
    pub mask_lowercase: bool,
    pub ranges: &'a ReportRanges<'a>,
    /// The run's program type (`RunRecord::program_type`, `in_pssm::program_type`): PSI-TBLASTN's
    /// results carry no query masks (`report::record::carries_query_masks`).
    pub program_type: Option<crate::report::record::ProgramType>,
}

/// Calls `each` with the queries and the archive of each batch of `spans`, in order
/// (`BatchArchive`: its request, `RunRequests::with_batch`, and its results).
fn for_each_batch(
    data: &ArchiveData<'_>,
    spans: &[Range<usize>],
    mut each: impl FnMut(Range<usize>, &BatchArchive<'_>) -> Result<()>,
) -> Result<()> {
    let requests = requests(
        matches!(data.shown, ShownSubjects::Fasta { .. }).then_some(data.subject_records),
        data.mask_lowercase,
    );
    for span in spans.iter().cloned() {
        if span.is_empty() {
            continue;
        }
        let records = &data.queries[span.clone()];
        let first = span.start;
        let searched = !data.query_batch_skipped[first];
        // Only a searched batch has results masks (the search made them): they are computed per
        // batch; a PSSM search's results have none.
        let masks = if searched && crate::report::record::carries_query_masks(data.program_type) {
            query_record_masks(records, first, data.seg, data.mask_lowercase, data.ranges)
        } else {
            Vec::new()
        };

        // The results' masks are those of the search, so only a searched batch has them: one
        // `notset` mask per query with its SEG and lower-case intervals (a protein query).
        let mut results_masks: Vec<MaskEntry> = Vec::new();
        if searched {
            for (offset, intervals) in masks.into_iter().enumerate() {
                if !intervals.is_empty() {
                    results_masks.push(MaskEntry {
                        id: records[offset].local_id.clone(),
                        best: records[offset].believed().map(|ids| ids.best.clone()),
                        frame: Frame::NotSet,
                        intervals,
                    });
                }
            }
        }

        let length = &data.parameters.lengths[first];
        let context = if data.query_validity[first] {
            FirstContext {
                karlin: Some((data.ungapped_karlin[first], data.gapped)),
                search_space: length.eff_searchsp,
                length_adjustment: length.length_adjustment,
            }
        } else {
            FirstContext {
                karlin: None,
                search_space: 0,
                length_adjustment: 0,
            }
        };
        let ka_blocks = match &context.karlin {
            Some((ungapped, gapped)) => [
                ka_block(ungapped, false),
                gapped
                    .as_ref()
                    .map_or(KaBlock::missing(true), |gapped| ka_block(gapped, true)),
            ],
            None => [KaBlock::missing(false), KaBlock::missing(true)],
        };
        requests.with_batch(records, data.ranges.queries.get(first), |request| {
            let batch = BatchArchive {
                request: *request,
                masks: &results_masks,
                ka_blocks,
                search_space: context.search_space,
                length_adjustment: context.length_adjustment,
                searched_space: batch_searched_space(
                    searched,
                    length,
                    &data.queries[first],
                    data.query_validity[first],
                ),
                round: None,
            };
            each(span, &batch)
        })?;
    }
    Ok(())
}

/// Writes the archive of each batch of `render` to `writer`; `warnings` writes the warnings of a
/// batch's queries (to the diagnostics, which saves them as the batch's messages) before its
/// archive, and its sink gets NCBI's posts of the bytes that the archive's stream replaces
/// (`write_archive`).
pub(super) fn write_tblastn_archives(
    handle: &ArchiveHandle,
    render: &ArchiveRender<'_>,
    data: &ArchiveData<'_>,
    warnings: &mut QueryWarnings<'_>,
    writer: &mut dyn crate::report::ostream::ReportWrite,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Tblastn,
        db_title: data.shown.db_title(),
        // NCBI's `CBlastFormatUtil`/`CBlastSeqalign` builds Std-seg diagonals for an ungapped
        // translated search (`report::serial::seqalign`).
        gapped: data.gapped.is_some(),
        program_type: None,
    };
    for_each_batch(data, render.spans, |span, batch| {
        for index in span.clone() {
            warnings.before_query(index, &mut *writer)?;
        }
        let first = span.start;
        let mut aligns = TblastnAligns {
            run: RunRecord {
                program: run.program,
                db_title: run.db_title,
                gapped: run.gapped,
                program_type: run.program_type,
            },
            queries: data.queries,
            first,
            hitlists: &data.hitlists[span.clone()],
            subject_records: data.subject_records,
            subject_letters: data.subject_letters,
            shown: data.shown,
            ranges: data.ranges,
        };
        write_archive(
            handle,
            batch,
            &mut aligns,
            &mut *writer,
            &mut *warnings.sink,
        )
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:123-128
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3800-3837
/// tblastn's part of `-import_search_strategy` ([`import::apply`]): the arguments its command
/// line keeps and its table of the strategy's setters over its own options.
impl import::StrategyArgs for TblastnArgs {
    const PROGRAM: Program = Program::Tblastn;

    fn runs(program: EProgram) -> bool {
        program == EProgram::Tblastn
    }

    fn imports(&self) -> bool {
        self.import_search_strategy.is_some() && self.import.is_none()
    }

    fn defaults() -> Result<Self> {
        match import::defaults_of("tblastn")? {
            crate::cli::Commands::Tblastn(args) => Ok(args),
            _ => unreachable!("the arguments of tblastn"),
        }
    }

    fn kept(self) -> import::Kept {
        let given = self.given.list().to_vec();
        import::Kept {
            path: self.import_search_strategy.expect("an import"),
            given,
            query: self.query,
            query_loc: self.query_loc,
            strand: None,
            lcase_masking: self.lcase_masking,
            parse_deflines: self.parse_deflines,
            id_scope: self.id_scope,
            out: self.out,
            outfmt: self.outfmt,
            show_gis: self.show_gis,
            html: self.html,
            num_descriptions: self.num_descriptions,
            num_alignments: self.num_alignments,
            max_target_seqs: self.max_target_seqs.map(|n| n as i64),
            formatting: self.formatting,
            num_threads: self.num_threads as i64,
            mt_mode: self.mt_mode,
            database: import::KeptDatabase::given(
                self.db,
                self.subject,
                self.subject_loc,
                self.db_filter,
                self.db_mask,
                self.dbsize,
            ),
            use_index: false,
            mb_index: import::MbIndexRequest::None,
        }
    }

    fn keep(&mut self, kept: import::Kept) -> Result<()> {
        self.query = kept.query;
        self.query_loc = kept.query_loc;
        self.lcase_masking = kept.lcase_masking;
        self.parse_deflines = kept.parse_deflines;
        self.id_scope = kept.id_scope;
        self.out = kept.out;
        self.outfmt = kept.outfmt;
        self.show_gis = kept.show_gis;
        self.html = kept.html;
        self.num_descriptions = kept.num_descriptions;
        self.num_alignments = kept.num_alignments;
        self.max_target_seqs = kept
            .max_target_seqs
            .map(|n| usize::try_from(n).unwrap_or(0));
        self.formatting = kept.formatting;
        self.num_threads = usize::try_from(kept.num_threads).unwrap_or(1);
        self.mt_mode = kept.mt_mode;
        Ok(())
    }

    fn fields(&mut self) -> import::Fields<'_> {
        import::Fields {
            import_search_strategy: &mut self.import_search_strategy,
            db: &mut self.db,
            subject: &mut self.subject,
            subject_loc: &mut self.subject_loc,
            db_filter: &mut self.db_filter,
            db_mask: &mut self.db_mask,
            dbsize: &mut self.dbsize,
            searchsp: &mut self.searchsp,
            formatting: &mut self.formatting,
            num_descriptions: self.num_descriptions,
            num_alignments: self.num_alignments,
            outfmt: &self.outfmt,
            import: &mut self.import,
        }
    }

    fn set_task(&mut self, task: &str) {
        self.task = task.to_string();
    }

    fn set(&mut self, setter: &Setter) -> Result<bool> {
        let program = Program::Tblastn.name();
        match setter {
            Setter::EvalueThreshold(v) => self.evalue = Some(*v),
            Setter::WordSize(n) => self.word_size = Some(*n),
            Setter::WordThreshold(n) => self.threshold = Some(f64::from(*n)),
            Setter::WindowSize(n) => self.window_size = Some(*n),
            Setter::MatrixName(name) => self.matrix = Some(name.clone().into()),
            Setter::GapOpeningCost(n) => self.gap_open = Some(*n),
            Setter::GapExtensionCost(n) => self.gap_extend = Some(*n),
            Setter::GapXDropoff(v) => self.xdrop_gap = Some(*v),
            Setter::GapXDropoffFinal(v) => self.xdrop_gap_final = Some(*v),
            Setter::CullingLimit(n) => self.culling_limit = Some(*n),
            Setter::MaxHspsPerSubject(n) => {
                self.max_hsps = usize::try_from(*n).ok().filter(|n| *n > 0)
            }
            Setter::QueryCovHspPerc(v) => self.qcov_hsp_perc = Some(*v),
            Setter::MaskAtHash(on) => self.soft_masking = Some(*on),
            Setter::GappedMode(gapped) => self.ungapped = !*gapped,
            Setter::SumStatisticsMode(on) => self.sum_stats = Some(*on),
            Setter::DbGeneticCode(n) => {
                self.db_gencode = import::Start::genetic_code(*n, "database", program)?
            }
            Setter::LongestIntronLength(n) => self.max_intron_length = Some(*n),
            Setter::SubjectBestHit => self.subject_besthit = true,
            Setter::BestHitOverhang(v) => self.best_hit_overhang = Some(*v),
            Setter::BestHitScoreEdge(v) => self.best_hit_score_edge = Some(*v),
            Setter::SmithWatermanMode(false) | Setter::OutOfFrameMode(false) => {}
            _ => return Ok(false),
        }
        Ok(true)
    }

    fn finish(&mut self, collected: import::Collected, _run: &import::ImportRun) -> Result<()> {
        let program = Program::Tblastn.name();
        collected.no_dust(program)?;
        collected.no_repeat(program)?;
        collected.no_nucleotide_algorithms(program)?;
        if let Some(value) = collected.seg.value() {
            self.seg = Some(value.into());
        }
        if let Some(value) = collected.composition_value(false, program)? {
            self.comp_based_stats = value.into();
        }
        Ok(())
    }

    fn into_command(self) -> crate::cli::Commands {
        crate::cli::Commands::Tblastn(self)
    }
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod tests;
