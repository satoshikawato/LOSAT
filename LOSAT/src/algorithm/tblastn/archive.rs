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

use std::io::Write;
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
use crate::report::archive::requests::{LcaseFrames, RequestPolicy, RunRequests};
use crate::report::archive::{
    ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive, KaBlock,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::serial::seqalign::{write_seq_aligns, Segments};
use crate::stats::tables::KarlinParams;
use crate::utils::seg::SegParams;

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
    inputs.matrix = args.matrix.clone();
    inputs.xdrop_gap = args.xdrop_gap;
    inputs.xdrop_gap_final = args.xdrop_gap_final;
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps.map(|count| count as i32);
    inputs.culling_limit = args.culling_limit;
    inputs.subject_besthit = args.subject_besthit;
    inputs.sum_stats = args.sum_stats;
    // `-threshold` is read as a double and recorded as an `int` (`static_cast<int>`), converted
    // as the x86-64 build does (`INT_MIN` past `Int4`, as the lookup table takes it:
    // `ncbi_int4_from_double`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:884-892
    inputs.threshold = args.threshold.map(ncbi_int4_from_double);
    inputs.window_size = args.window_size;
    inputs.ungapped = args.ungapped;
    inputs.soft_masking = args.soft_masking;
    inputs.seg = args
        .seg
        .as_deref()
        .and_then(crate::report::archive::options::seg_filter);
    inputs.max_intron_length = args.max_intron_length.unwrap_or(0);
    inputs.db_gencode = i32::from(args.db_gencode);
    inputs.imported = args.import.as_ref().map(|run| run.recorded.clone());
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(args.formatting.recorded_hitlist_start(), |size| size as i32);
    if let Ok((mode, _)) =
        crate::blastinput::app::parse_comp_based_stats(&args.comp_based_stats, false, false)
    {
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
        args.database_name(),
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
    /// (both blocks are then missing).
    karlin: Option<(KarlinParams, KarlinParams)>,
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
        out: &mut ArchiveOut<'_, E>,
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
            out.write_with(|enc, flush| {
                write_seq_aligns(
                    enc,
                    &self.run,
                    &query,
                    &subjects,
                    &hsps,
                    &mut segments,
                    |enc| flush(enc),
                );
            });
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
    pub gapped: KarlinParams,
    pub seg: Option<&'a SegParams>,
    pub mask_lowercase: bool,
    pub ranges: &'a ReportRanges<'a>,
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
        // batch.
        let masks = if searched {
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
            Some((ungapped, gapped)) => [ka_block(ungapped, false), ka_block(gapped, true)],
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
    writer: &mut dyn Write,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Tblastn,
        db_title: data.shown.db_title(),
        gapped: true,
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

/// `-import_search_strategy` of TBLASTN (`RecoverSearchStrategy`, `SetOptionsForSavedStrategy`,
/// blast_args.cpp:3800-3836, quoted at `blastp::archive::apply_strategy`): the arguments of a
/// command line that gives nothing (the strategy task's defaults) with the strategy's options set
/// on them in the file's order, and from `args` only what NCBI's re-read argument classes take:
/// the query options, `-out`, the formatting options, the threads and the database arguments
/// when `-db` or `-subject` is given. The run then goes on as any other (one search path).
pub(crate) fn apply_strategy(
    args: super::args::TblastnArgs,
    diagnostics: &mut dyn std::io::Write,
) -> Result<super::args::TblastnArgs> {
    use crate::blastinput::app::Program;
    use crate::report::archive::read::import::{self, Command, Start};
    use crate::report::archive::read::options::{EProgram, Setter};
    let program = Program::Tblastn.name();
    let path = args.import_search_strategy.clone().expect("an import");
    let source = args.import.as_ref().and_then(|run| run.source.clone());
    let start = import::start(
        &Command {
            path: &path,
            bytes: source.as_deref().map(Vec::as_slice),
            program,
            override_query: args.query.as_os_str() != "-",
            override_subject: args.db.is_some() || args.subject.is_some(),
            use_index_set: false,
        },
        "tblastn",
        diagnostics,
    )?;
    if start.imported.program() != EProgram::Tblastn {
        return Err(start.unported_program(program));
    }
    let cli: crate::cli::Cli = crate::cli::try_parse_from(["losat", "tblastn"])?;
    let crate::cli::Commands::Tblastn(mut out) = cli.command else {
        unreachable!("the arguments of tblastn")
    };
    let mut run = start.run(args.parse_deflines, program)?;
    run.query_range = start.query_range(args.query_loc.as_ref(), program)?;
    // `CQueryOptionsArgs`, `CStdCmdLineArgs`, `CFormattingArgs`, `CMTArgs`.
    out.query = args.query;
    out.query_loc = args.query_loc;
    out.lcase_masking = args.lcase_masking;
    out.parse_deflines = args.parse_deflines;
    out.id_scope = args.id_scope;
    out.out = args.out;
    out.outfmt = args.outfmt;
    out.show_gis = args.show_gis;
    out.num_descriptions = args.num_descriptions;
    out.num_alignments = args.num_alignments;
    out.max_target_seqs = args.max_target_seqs;
    out.formatting = args.formatting;
    out.num_threads = args.num_threads;
    out.mt_mode = args.mt_mode;
    out.import_search_strategy = Some(path);
    // `CBlastDatabaseArgs` when the command line gives the database or the subjects, else the
    // strategy's.
    // `-dbsize` with `-db` or `-subject` is set after the strategy's options, unless the
    // strategy's search space trumps it (blast_args.cpp:2565-2572).
    let mut cli_dbsize = None;
    if args.db.is_some() || args.subject.is_some() {
        out.db = args.db;
        out.subject = args.subject;
        out.subject_loc = args.subject_loc;
        out.db_filter = args.db_filter;
        cli_dbsize = args.dbsize;
    } else if let Some(name) = start.database() {
        out.db = Some(name.to_string());
        out.db_filter = start.db_filter(program)?;
        if start.db_mask().is_some() {
            return Err(Start::unported(
                "a search strategy with a database mask",
                program,
            ));
        }
    } else if start.has_subjects() {
        // The Owner's decision F2-Q1 (2026-11-20 「Q1=C」): a strategy's subjects are searched
        // once, as with `-num_threads 1`, whatever `-num_threads` says (NCBI searches them once per
        // thread: its thread check reads only the command line's `-subject`).
        out.num_threads = 1;
    }
    out.task = start.imported.program().task_name().to_string();
    let options = &start.imported.options.options;
    let mut seg = import::SegSetters::default();
    let (mut cbs, mut unified_p) = (None, 0);
    for setter in &options.setters {
        if seg.take(setter) {
            continue;
        }
        match setter {
            Setter::EvalueThreshold(v) => out.evalue = Some(*v),
            Setter::WordSize(n) => out.word_size = Some(*n),
            Setter::WordThreshold(n) => out.threshold = Some(f64::from(*n)),
            Setter::WindowSize(n) => out.window_size = Some(*n),
            Setter::MatrixName(name) => out.matrix = Some(name.clone()),
            Setter::GapOpeningCost(n) => out.gap_open = Some(*n),
            Setter::GapExtensionCost(n) => out.gap_extend = Some(*n),
            Setter::GapXDropoff(v) => out.xdrop_gap = Some(*v),
            Setter::GapXDropoffFinal(v) => out.xdrop_gap_final = Some(*v),
            Setter::CompositionBasedStats(n) => cbs = Some(*n),
            Setter::UnifiedP(n) => unified_p = *n,
            Setter::CullingLimit(n) => out.culling_limit = Some(*n),
            Setter::EffectiveSearchSpace(n) => {
                out.searchsp = Some(*n);
                run.searchsp_from_strategy = true;
            }
            Setter::DbLength(n) => out.dbsize = Some(*n),
            Setter::HitlistSize(n) => out.formatting.hitlist_start = usize::try_from(*n).ok(),
            Setter::MaxHspsPerSubject(n) => {
                out.max_hsps = usize::try_from(*n).ok().filter(|n| *n > 0)
            }
            Setter::QueryCovHspPerc(v) => out.qcov_hsp_perc = Some(*v),
            Setter::MaskAtHash(on) => out.soft_masking = Some(*on),
            Setter::GappedMode(gapped) => out.ungapped = !*gapped,
            Setter::SumStatisticsMode(on) => out.sum_stats = Some(*on),
            Setter::DbGeneticCode(n) => {
                out.db_gencode = u8::try_from(*n).map_err(|_| {
                    Start::unported("a search strategy with this database genetic code", program)
                })?
            }
            Setter::LongestIntronLength(n) => out.max_intron_length = Some(*n),
            Setter::DustFiltering(false)
            | Setter::RepeatFiltering(false)
            | Setter::SmithWatermanMode(false)
            | Setter::OutOfFrameMode(false) => {}
            other => return Err(Start::unported_setter(other, program)),
        }
    }
    if let Some(value) = seg.value() {
        out.seg = Some(value);
    }
    if let Some(mode) = cbs {
        if !(0..=3).contains(&mode) || unified_p != 0 {
            return Err(Start::unported(
                &format!("a search strategy with composition-based statistics {mode} (unified P {unified_p})"),
                program,
            ));
        }
        out.comp_based_stats = mode.to_string();
    }
    if out.searchsp.is_none() {
        if let Some(size) = cli_dbsize {
            out.dbsize = Some(size);
        }
    }
    import::formatting_early_return(
        &mut out.formatting,
        out.num_descriptions.is_some(),
        out.num_alignments.is_some(),
        &out.outfmt,
    );
    out.import = Some(run);
    Ok(out)
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod tests;
