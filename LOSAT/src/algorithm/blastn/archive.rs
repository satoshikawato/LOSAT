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
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(500, |size| size as i32);
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
                    round: None,
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
        pssm_type: None,
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

#[cfg(test)]
#[path = "archive_tests.rs"]
mod tests;
