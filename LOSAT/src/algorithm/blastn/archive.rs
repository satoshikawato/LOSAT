//! BLASTN's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:302-307
//! One archive per query batch, written after the batch's search; `report::archive` writes
//! the document, this module gathers what it holds from the hit lists, the query masks and
//! the Karlin-Altschul blocks that BLASTN keeps for its other formats.

use std::io::Write;
use std::ops::Range;

use anyhow::Result;

use super::args::BlastnArgs;
use super::coordination::collect_lowercase_masks;
use super::hsp::BlastnHitList;
use super::pairwise::shown_query_masks;
use super::scoring::ContextKarlin;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::blastinput::value_parsers::{parse_dust_filtering, DustSpec};
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{Filter, OptionInputs};
use crate::report::archive::requests::{LcaseFrames, RequestPolicy, RunRequests};
use crate::report::archive::{
    ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive, KaBlock,
    RequestInputs,
};
#[cfg(test)]
use crate::report::diag_log::ERROR;
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryNames, RunRecord, SubjectNames, SubjectRecord};
use crate::report::serial::seqalign::{write_seq_aligns, Segments};
use crate::utils::dust::MaskedInterval;

/// The recorded options of a BLASTN run (`options::algorithm_options`'s inputs). `dust` is
/// the `-dust` value as given (`BlastnArgs::resolve_dust` consumes it).
pub(crate) fn option_inputs(args: &BlastnArgs, dust: Option<&str>) -> OptionInputs {
    let mut inputs = OptionInputs::new(Program::Blastn, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    inputs.evalue = args.evalue;
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.word_size = args.word_size.map(|size| size as i32);
    inputs.perc_identity = args.percent_identity;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    // (An ungapped search's archive: `UngappedMode TRUE`, no gapped `ka-block`, Dense-diags;
    // the search strategy records it too, probe `ung` of RG-F1.)
    inputs.ungapped = args.ungapped;
    inputs.max_hsps = args.max_hsps_per_subject.map(|count| count as i32);
    inputs.penalty = args.penalty;
    inputs.reward = args.reward;
    inputs.template_type = args.template_type.map(|kind| {
        use crate::algorithm::blastn::disc_lookup::DiscWordType;
        match kind {
            DiscWordType::Coding => 0,
            DiscWordType::Optimal => 1,
            DiscWordType::TwoTemplates => 2,
        }
    });
    inputs.template_length = args.template_length.map(i32::from);
    inputs.dust = dust
        .and_then(|value| parse_dust_filtering(value).ok())
        .map(|spec| match spec {
            DustSpec::No => Filter::No,
            DustSpec::Yes => Filter::Yes,
            DustSpec::Parameters {
                level,
                window,
                linker,
            } => Filter::Dust {
                level: level as i32,
                window: window as u32 as i32,
                linker: linker as u32 as i32,
            },
        });
    inputs.soft_masking = Some(args.soft_masking());
    inputs.filtering_db = args.filtering_db.clone();
    inputs.window_masker_db = args.window_masker_db.clone();
    inputs.window_masker_taxid = args.window_masker_taxid;
    inputs.culling_limit = args.culling_limit;
    inputs.best_hit_overhang = args.best_hit_overhang;
    inputs.best_hit_score_edge = args.best_hit_score_edge;
    inputs.subject_besthit = args.subject_besthit;
    // The X-drops and the window as given (CGenericSearchArgs and CWindowSizeArg set them
    // only then, blast_args.cpp:280-286,487-489) and `-no_greedy` (the algorithms); the
    // ungapped X-drop and the off-diagonal range are not recorded (`options.rs`).
    inputs.xdrop_gap = args.xdrop_gap;
    inputs.xdrop_gap_final = args.xdrop_gap_final;
    inputs.window_size = args.window_size;
    inputs.no_greedy = args.no_greedy;
    inputs.mb_index_name = args
        .megablast_index_name()
        .map(|(name, _)| name.to_string());
    inputs.imported = args.import.as_ref().map(|run| run.recorded.clone());
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(args.formatting.recorded_hitlist_start(), |size| size as i32);
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with BLASTN's recorded options).
pub(crate) fn handle_for(
    args: &BlastnArgs,
    dust: Option<&str>,
    archive_requested: bool,
    options: &crate::report::archive::OptionsHandle,
) -> Option<ArchiveHandle> {
    if !archive_requested {
        return None;
    }
    Some(ArchiveHandle::for_program(
        &option_inputs(args, dust),
        args.database_name(),
        options,
    ))
}

/// The requests of a run's batches (`RunRequests`): nucleotide queries and subjects (the
/// `-subject` records `subjects`, `None` for `-db`), and with `-lcase_masking` (`lcase_masking`)
/// the lower-case runs of each query for both strands. Before the search (the search strategy's
/// requests) the query vector's mask has no DUST or WindowMasker interval yet, which
/// `x_CalculateMasks` adds in the search (probe `lcq` of RG-F1); the archive writes its own
/// `LCaseMask` parameters (`for_each_batch`).
pub(crate) fn requests(subjects: Option<&[FastaRecord]>, lcase_masking: bool) -> RunRequests<'_> {
    RunRequests::new(
        RequestPolicy {
            query_protein: false,
            subject_protein: false,
            lcase: lcase_masking.then_some(LcaseFrames::PlusTwice),
        },
        subjects,
    )
}

pub(crate) use crate::report::archive::error_lines;

/// The queries of a call's batches and what their archives read (`for_each_batch`).
pub(crate) struct CallBatches<'a> {
    /// The call's queries: one batch, or several when the reports of every batch are written
    /// after the search.
    pub queries: &'a [FastaRecord],
    pub query_lengths: &'a [usize],
    pub query_placements: &'a Placements,
    /// `-subject`: every subject record, whole (the Bioseqs of the request).
    pub subject_records: Option<&'a [FastaRecord]>,
    pub query_masks: &'a [Vec<MaskedInterval>],
    /// The DUST-only masks of the call's queries (plus strand, in the searched letters),
    /// made for every query; with `-window_masker_db` (`window_masker`) the WindowMasker masks.
    /// Only the archive reads them (empty without it).
    pub dust_masks: &'a [Vec<MaskedInterval>],
    /// With `-filtering_db`, the repeat hits of the call's queries (0-based, inclusive, in the
    /// searched letters; `repeat_filter::RepeatFilter::hits`), empty for a query without one.
    /// Only the archive reads them (empty without it or without the option).
    pub repeat_hits: &'a [Vec<(usize, usize)>],
    pub unsearched: &'a [bool],
    pub lcase_masking: bool,
    pub window_masker: bool,
    pub both_strands: bool,
    pub minus_strand: bool,
    /// The queries of each batch of the call, in order; empty for a call that is one batch.
    pub spans: &'a [Range<usize>],
    /// The search is not `-ungapped`: an ungapped search has no gapped Karlin-Altschul block
    /// (`ka-blocks`) and its Seq-aligns are Dense-diags (`RunRecord::gapped`).
    pub gapped: bool,
}

/// What the archive of a call needs: the call's batches and their results.
pub(crate) struct ArchiveCall<'a> {
    pub handle: &'a ArchiveHandle,
    pub batches: CallBatches<'a>,
    /// The subject names of the Seq-aligns.
    pub subjects: SubjectNames<'a>,
    pub subject_placements: &'a Placements,
    pub hit_lists: &'a [Option<BlastnHitList>],
}

/// The values of a batch's first query that `search-stats`, `ka-blocks` and the
/// `EffectiveSearchSpace` option take.
pub(crate) struct FirstContext {
    /// The ungapped and gapped blocks of the first valid context; `None` when the first
    /// query has no valid context.
    pub karlin: Option<ContextKarlin>,
    /// The effective search space and length adjustment of that context (0 without one).
    pub search_space: i64,
    pub length_adjustment: i64,
    /// The search space of the batch's context 0 (the option).
    pub context0_space: i64,
    /// That context's query length, 0 for an invalid one (`archive::RequestLength::space`).
    pub context0_length: i64,
}

/// The Seq-aligns of one batch of a call.
struct BlastnAligns<'a> {
    run: RunRecord<'a>,
    queries: QueryNames<'a>,
    subjects: SubjectNames<'a>,
    subject_placements: &'a Placements,
    hit_lists: &'a [Option<BlastnHitList>],
    /// The index of the batch's first query in the call.
    first: usize,
}

impl<'a> AlignSource for BlastnAligns<'a> {
    fn write<E: crate::report::serial::encoder::Encoder>(
        &mut self,
        out: &mut ArchiveOut<'_, E>,
    ) -> Result<()> {
        let mut segments = Segments::default();
        let mut subject_records: Vec<SubjectRecord<'a>> = Vec::new();
        let mut hsp_records: Vec<HspRecord<'a>> = Vec::new();
        let hit_lists: &'a [Option<BlastnHitList>] = self.hit_lists;
        for (offset, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            super::seqalign::query_records(
                hit_list,
                &self.subjects,
                self.subject_placements,
                &mut subject_records,
                &mut hsp_records,
            );
            let query = self.queries.record(self.first + offset);
            out.write_with(|enc, flush| {
                write_seq_aligns(
                    enc,
                    &self.run,
                    &query,
                    &subject_records,
                    &hsp_records,
                    &mut segments,
                    |enc| flush(enc),
                );
            });
        }
        Ok(())
    }
}

/// The intervals of a query's `LCaseMask` parameter, in record coordinates.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:110-127
/// The query vector's mask is the lower case of the whole record (`-lcase_masking`, read
/// before any `-query_loc` cut, for both strands) merged with the DUST intervals of the
/// searched range, which the mapper puts in the order of the query's strand. A merge
/// sorts and joins overlapping or abutting intervals; a DUST-only mask is merged
/// (`Merge`) in the order the mapper gave: descending for `-strand minus`. Each interval
/// is written once for each strand of the mask (twice with a lower-case mask or both
/// strands).
///
/// With `-window_masker_db` (`window_masker`) the query vector's mask is WindowMasker's,
/// which replaced the others: every interval once for each frame, in the masker's order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:179-189
fn lcase_mask_intervals(
    record: &FastaRecord,
    dust: &[MaskedInterval],
    offset: usize,
    lcase_masking: bool,
    window_masker: bool,
    both_strands: bool,
    minus_strand: bool,
) -> Vec<(i32, i32)> {
    if window_masker {
        return dust
            .iter()
            .flat_map(|mask| {
                let interval = ((mask.start + offset) as i32, (mask.end + offset) as i32 - 1);
                [interval, interval]
            })
            .collect();
    }
    let user: Vec<(i32, i32)> = if lcase_masking {
        collect_lowercase_masks(record.archive_letters())
            .iter()
            .map(|mask| (mask.start as i32, mask.end as i32 - 1))
            .collect()
    } else {
        Vec::new()
    };
    let user_mask = !user.is_empty();
    let mut all = user;
    all.extend(
        dust.iter()
            .map(|mask| ((mask.start + offset) as i32, (mask.end + offset) as i32 - 1)),
    );
    all.sort_unstable();
    let mut merged: Vec<(i32, i32)> = Vec::with_capacity(all.len());
    for (from, to) in all {
        match merged.last_mut() {
            Some(last) if from <= last.1 + 1 => last.1 = last.1.max(to),
            _ => merged.push((from, to)),
        }
    }
    if minus_strand && !user_mask {
        merged.reverse();
    }
    let times = if both_strands || user_mask { 2 } else { 1 };
    merged
        .iter()
        .flat_map(|&interval| std::iter::repeat_n(interval, times))
        .collect()
}

/// The `LCaseMask` intervals of a query whose repeat search (`-filtering_db`) found `hits` (in the
/// searched letters, which start at `offset`): the query vector's masks before the search (the
/// lower case of the whole record with `-lcase_masking`, DUST's of the searched letters) combined
/// with the hits (`repeat_filter::combine`), in ascending order, each interval for both strands,
/// whatever `-strand` is.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:290-304
fn repeat_mask_intervals(
    record: &FastaRecord,
    dust: &[MaskedInterval],
    offset: usize,
    lcase_masking: bool,
    hits: &[(usize, usize)],
) -> Vec<(i32, i32)> {
    let mut all: Vec<(usize, usize)> = if lcase_masking {
        collect_lowercase_masks(record.archive_letters())
            .iter()
            .map(|mask| (mask.start, mask.end - 1))
            .collect()
    } else {
        Vec::new()
    };
    all.extend(
        dust.iter()
            .map(|mask| (mask.start + offset, mask.end + offset - 1)),
    );
    all.extend(hits.iter().map(|&(from, to)| (from + offset, to + offset)));
    super::repeat_filter::combine(all)
        .into_iter()
        .flat_map(|(from, to)| {
            let interval = (from as i32, to as i32);
            [interval, interval]
        })
        .collect()
}

/// Calls `each` with the queries and the archive of each batch of `call`, in order
/// (`BatchArchive`: its request, `RunRequests::with_batch`, and its results).
/// `first_context(first_query)` gives the values of a batch's first query.
///
/// The `LCaseMask` parameters are the masks of NCBI's query vector after the search, which hold
/// the DUST or WindowMasker intervals that the search added (`lcase_mask_intervals`), so the
/// archive writes its own (the search strategy's requests, written before the search, hold the
/// lower-case runs only, `requests`).
pub(crate) fn for_each_batch(
    call: &CallBatches<'_>,
    first_context: &dyn Fn(usize) -> Result<FirstContext>,
    mut each: impl FnMut(Range<usize>, &BatchArchive<'_>) -> Result<()>,
) -> Result<()> {
    let whole = [0..call.queries.len()];
    let spans: &[Range<usize>] = if call.spans.is_empty() {
        &whole
    } else {
        call.spans
    };
    // The request without `LCaseMask` parameters: the archive writes its own (above).
    let requests = requests(call.subject_records, false);
    for span in spans.iter().cloned() {
        if span.is_empty() {
            continue;
        }
        let first = span.start;
        let searched = !call.unsearched[first];
        let context = first_context(first)?;

        let mut masks: Vec<MaskEntry> = Vec::new();
        let mut lcase_masks: Vec<MaskEntry> = Vec::new();
        for index in span.clone() {
            let record = &call.queries[index];
            let offset = call.query_placements.offset(index);
            // The results' masks are those of the search, so only a searched batch has them
            // (and leave out a mask over the whole query, `shown_query_masks`).
            if searched {
                let shown: Vec<(i32, i32)> =
                    shown_query_masks(&call.query_masks[index], call.query_lengths[index], offset)
                        .iter()
                        .map(|mask| (mask.start as i32, mask.end as i32 - 1))
                        .collect();
                if !shown.is_empty() {
                    masks.push(MaskEntry {
                        id: record.local_id.clone(),
                        best: record.believed().map(|ids| ids.best.clone()),
                        frame: Frame::Plus1,
                        intervals: shown,
                    });
                }
            }
            // The `LCaseMask` parameters are the masks of NCBI's query vector, which every
            // query has, valid or not and searched or not (`x_CalculateMasks`).
            let intervals = match call.repeat_hits.get(index).filter(|hits| !hits.is_empty()) {
                Some(hits) => repeat_mask_intervals(
                    record,
                    &call.dust_masks[index],
                    offset,
                    call.lcase_masking,
                    hits,
                ),
                None => lcase_mask_intervals(
                    record,
                    &call.dust_masks[index],
                    offset,
                    call.lcase_masking,
                    call.window_masker,
                    call.both_strands,
                    call.minus_strand,
                ),
            };
            if !intervals.is_empty() {
                lcase_masks.push(MaskEntry {
                    id: record.local_id.clone(),
                    best: record.believed().map(|ids| ids.best.clone()),
                    frame: Frame::Plus1,
                    intervals,
                });
            }
        }

        let ka_blocks = match &context.karlin {
            Some(blocks) => [
                ka_block(&blocks.ungapped, false),
                // An ungapped search has no gapped block (`KaBlock::missing`, as BLASTX's).
                if call.gapped {
                    ka_block(&blocks.gapped, true)
                } else {
                    KaBlock::missing(true)
                },
            ],
            None => [KaBlock::missing(false), KaBlock::missing(true)],
        };
        requests.with_batch(
            &call.queries[span.clone()],
            call.query_placements.get(first),
            |request| {
                let batch = BatchArchive {
                    request: RequestInputs {
                        lcase_masks: &lcase_masks,
                        ..*request
                    },
                    masks: &masks,
                    ka_blocks,
                    search_space: context.search_space,
                    length_adjustment: context.length_adjustment,
                    searched_space: searched.then_some(crate::report::archive::SearchedContext0 {
                        space: context.context0_space,
                        query_length: context.context0_length,
                    }),
                };
                each(span, &batch)
            },
        )?;
    }
    Ok(())
}

/// Writes the archive of each batch of `call` to `writer`; `warnings` writes the warnings
/// of a batch's queries (to the diagnostics, which saves them as the batch's messages)
/// before its archive, and its sink gets NCBI's posts of the bytes that the archive's stream
/// replaces (`write_archive`). `first_context(first_query)` gives the values of a batch's first
/// query.
pub(crate) fn write_blastn_archives(
    call: &ArchiveCall<'_>,
    warnings: &mut QueryWarnings<'_>,
    first_context: &dyn Fn(usize) -> Result<FirstContext>,
    writer: &mut dyn Write,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Blastn,
        db_title: call.subjects.db_title(),
        gapped: call.batches.gapped,
    };
    for_each_batch(&call.batches, first_context, |span, batch| {
        for index in span.clone() {
            warnings.before_query(index, &mut *writer)?;
        }
        let mut aligns = BlastnAligns {
            run,
            queries: QueryNames {
                records: call.batches.queries,
                searched_lengths: call.batches.query_lengths,
                placements: call.batches.query_placements,
            },
            subjects: call.subjects,
            subject_placements: call.subject_placements,
            hit_lists: &call.hit_lists[span.clone()],
            first: span.start,
        };
        write_archive(
            call.handle,
            batch,
            &mut aligns,
            &mut *writer,
            &mut *warnings.sink,
        )
    })
}

/// The arguments a blastn run keeps for another program's strategy (`import::RunAs`).
fn carried(args: BlastnArgs) -> crate::report::archive::read::import::Carried {
    crate::report::archive::read::import::Carried {
        app: "blastn",
        path: args.import_search_strategy.expect("an import"),
        source: args.import.as_ref().and_then(|run| run.source.clone()),
        query: args.query,
        query_loc: args.query_loc,
        strand: args.strand,
        lcase_masking: args.lcase_masking,
        parse_deflines: args.parse_deflines,
        id_scope: args.id_scope,
        out: args.out,
        outfmt: args.outfmt,
        show_gis: args.show_gis,
        num_descriptions: args.num_descriptions,
        num_alignments: args.num_alignments,
        max_target_seqs: args.max_target_seqs.map(|n| n as i64),
        formatting: args.formatting,
        num_threads: args.num_threads,
        mt_mode: args.mt_mode,
        db: args.db,
        subject: args.subject,
        subject_loc: args.subject_loc,
        db_filter: args.db_filter,
        dbsize: args.dbsize,
    }
}

/// blastn's arguments for a BLASTN strategy that another program runs (`import::RunAs`): its
/// defaults with the arguments the running program keeps.
pub(crate) fn from_carried(
    carried: crate::report::archive::read::import::Carried,
) -> Result<BlastnArgs> {
    let cli: crate::cli::Cli = crate::cli::try_parse_from(["losat", "blastn"])?;
    let crate::cli::Commands::Blastn(mut out) = cli.command else {
        unreachable!("the arguments of blastn")
    };
    out.import = carried.import_app();
    out.max_target_seqs = carried
        .max_target_seqs
        .map(|n| usize::try_from(n).unwrap_or(0));
    out.import_search_strategy = Some(carried.path);
    out.query = carried.query;
    out.query_loc = carried.query_loc;
    out.strand = carried.strand;
    out.lcase_masking = carried.lcase_masking;
    out.parse_deflines = carried.parse_deflines;
    out.id_scope = carried.id_scope;
    out.out = carried.out;
    out.outfmt = carried.outfmt;
    out.show_gis = carried.show_gis;
    out.num_descriptions = carried.num_descriptions;
    out.num_alignments = carried.num_alignments;
    out.formatting = carried.formatting;
    out.num_threads = carried.num_threads;
    out.mt_mode = carried.mt_mode;
    out.db = carried.db;
    out.subject = carried.subject;
    out.subject_loc = carried.subject_loc;
    out.db_filter = carried.db_filter;
    out.dbsize = carried.dbsize;
    out.command_line = Some(crate::cli::command_line_bytes());
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3810-3830
/// `-import_search_strategy` of BLASTN (`RecoverSearchStrategy`, `SetOptionsForSavedStrategy`):
/// the arguments of a command line that gives nothing (the strategy task's defaults) with the
/// strategy's options set on them in the file's order, and from `args` only what NCBI's
/// re-read argument classes take: the query options, `-out`, the formatting options, the
/// threads, the database arguments when `-db` or `-subject` is given and `-use_index` when
/// given. The run then goes on as any other (one search path).
pub(crate) fn apply_strategy(
    args: BlastnArgs,
    diagnostics: &mut dyn std::io::Write,
) -> Result<BlastnArgs> {
    use crate::report::archive::read::import::{self, Command, DbMaskAlgorithm, Start};
    use crate::report::archive::read::options::{EProgram, Setter};
    let program = Program::Blastn.name();
    let path = args.import_search_strategy.clone().expect("an import");
    let source = args.import.as_ref().and_then(|run| run.source.clone());
    let start = import::start(
        &Command {
            path: &path,
            bytes: source.as_deref().map(Vec::as_slice),
            program,
            override_query: args.query.as_os_str() != "-",
            override_subject: args.db.is_some() || args.subject.is_some(),
            use_index_set: args.use_index == Some(true),
        },
        args.import.as_ref().map_or("blastn", |run| run.app),
        diagnostics,
    )?;
    if !matches!(
        start.imported.program(),
        EProgram::Blastn | EProgram::Megablast | EProgram::DiscMegablast
    ) {
        return Err(start.run_as(program, || carried(args)));
    }
    let cli: crate::cli::Cli = crate::cli::try_parse_from(["losat", "blastn"])?;
    let crate::cli::Commands::Blastn(mut out) = cli.command else {
        unreachable!("the arguments of blastn")
    };
    // `CMbIndexArgs::ExtractAlgorithmOptions` (blast_args.cpp:3386-3412): `-use_index true` of
    // task megablast (the command line's `-task`) names its index by `-db`, which the strategy's
    // database is not (probe `I.blastn.db.useidx`).
    if args.use_index == Some(true) && args.task == "megablast" && args.db.is_none() {
        start.post_warnings();
        return Err(crate::blastinput::app::options_error(
            "Can not deduce database index name",
        ));
    }
    let mut run = start.run(args.parse_deflines, program)?;
    run.query_range = start.query_range(args.query_loc.as_ref(), program)?;
    // `CQueryOptionsArgs`, `CStdCmdLineArgs`, `CFormattingArgs`, `CMTArgs`, `CMbIndexArgs`.
    out.query = args.query;
    out.query_loc = args.query_loc;
    out.strand = args.strand;
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
    out.use_index = args.use_index;
    out.command_line = args.command_line;
    out.progress_bar = args.progress_bar;
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
        out.db_mask = args.db_mask;
        cli_dbsize = args.dbsize;
    } else if let Some(name) = start.database() {
        out.db = Some(name.to_string());
        out.db_filter = start.db_filter(program)?;
        if let Some((algorithm, kind)) = start.db_mask() {
            let value = match algorithm {
                DbMaskAlgorithm::Id(id) => id.to_string(),
                DbMaskAlgorithm::Key(key) => key,
            };
            match kind {
                1 => out.db_mask.db_soft_mask = Some(value),
                2 => out.db_mask.db_hard_mask = Some(value),
                _ => {
                    return Err(Start::unported(
                        "a search strategy with this database masking type",
                        program,
                    ))
                }
            }
        }
    } else if start.has_subjects() {
        // The Owner's decision F2-Q1 (2026-11-20 「Q1=C」): a strategy's subjects are searched
        // once, as with `-num_threads 1`, whatever `-num_threads` says (NCBI searches them once per
        // thread or aborts: its thread check reads only the command line's `-subject`).
        out.num_threads = 1;
    }
    out.task = start.imported.program().task_name().to_string();
    let options = &start.imported.options.options;
    let mut dust = import::DustSetters::default();
    let (mut extension, mut traceback) = (None, None);
    let (mut template_type, mut template_length) = (None, None);
    for setter in &options.setters {
        if dust.take(setter) {
            continue;
        }
        match setter {
            Setter::EvalueThreshold(v) => out.evalue = Some(*v),
            Setter::WordSize(n) => out.word_size = usize::try_from(*n).ok(),
            Setter::WindowSize(n) => out.window_size = Some(*n),
            Setter::MatchReward(n) => out.reward = Some(*n),
            Setter::MismatchPenalty(n) => out.penalty = Some(*n),
            Setter::GapOpeningCost(n) => out.gap_open = Some(*n),
            Setter::GapExtensionCost(n) => out.gap_extend = Some(*n),
            Setter::GapXDropoff(v) => out.xdrop_gap = Some(*v),
            Setter::GapXDropoffFinal(v) => out.xdrop_gap_final = Some(*v),
            Setter::GapExtnAlgorithm(n) => extension = Some(*n),
            Setter::GapTracebackAlgorithm(n) => traceback = Some(*n),
            Setter::CullingLimit(n) => out.culling_limit = Some(*n),
            Setter::EffectiveSearchSpace(n) => out.searchsp = Some(*n),
            Setter::DbLength(n) => out.dbsize = Some(*n),
            Setter::HitlistSize(n) => out.formatting.hitlist_start = usize::try_from(*n).ok(),
            Setter::MaxHspsPerSubject(n) => {
                out.max_hsps_per_subject = usize::try_from(*n).ok().filter(|n| *n > 0)
            }
            Setter::QueryCovHspPerc(v) => out.qcov_hsp_perc = Some(*v),
            Setter::PercentIdentity(v) => out.percent_identity = Some(*v),
            Setter::SubjectBestHit => out.subject_besthit = true,
            Setter::BestHitOverhang(v) => out.best_hit_overhang = Some(*v),
            Setter::BestHitScoreEdge(v) => out.best_hit_score_edge = Some(*v),
            Setter::MaskAtHash(on) => out.soft_masking = Some(*on),
            Setter::GappedMode(gapped) => out.ungapped = !*gapped,
            Setter::WindowMaskerTaxId(n) => out.window_masker_taxid = Some(*n),
            Setter::MbTemplateType(n) => template_type = Some(*n),
            Setter::MbTemplateLength(n) => template_length = Some(*n),
            Setter::SegFiltering(false) | Setter::RepeatFiltering(false) => {}
            other => return Err(Start::unported_setter(other, program)),
        }
    }
    if let Some(value) = dust.value() {
        out.dust_filtering = Some(value);
    }
    // `-no_greedy` sets both algorithms to dynamic programming (blast_args.cpp:701-704); the
    // greedy pair is megablast's default.
    match (extension, traceback) {
        (None, None) => {}
        (Some(0), Some(0)) => out.no_greedy = true,
        (Some(1), Some(1)) if out.task == "megablast" => {}
        _ => {
            return Err(Start::unported(
                "a search strategy with this pair of gapped extension algorithms",
                program,
            ))
        }
    }
    match (template_type, template_length) {
        (None, None) => {}
        (Some(kind), Some(length)) => {
            out.template_type = Some(match kind {
                0 => super::disc_lookup::DiscWordType::Coding,
                1 => super::disc_lookup::DiscWordType::Optimal,
                2 => super::disc_lookup::DiscWordType::TwoTemplates,
                _ => {
                    return Err(Start::unported(
                        "a search strategy with this discontiguous template type",
                        program,
                    ))
                }
            });
            out.template_length = Some(length);
        }
        _ => {
            return Err(Start::unported(
                "a search strategy with a template type or length alone",
                program,
            ))
        }
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
