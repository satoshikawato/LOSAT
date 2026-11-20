//! BLASTN's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:302-306
//! ```c++
//!             if (isArchiveFormat) {
//!                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//!                 m_Bah.ResetMessages();
//!             } else {
//! ```
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
    inputs.window_masker_db = args.window_masker_db.clone();
    inputs.window_masker_taxid = args.window_masker_taxid;
    inputs.culling_limit = args.culling_limit;
    inputs.subject_besthit = args.subject_besthit;
    // The X-drops and the window as given (CGenericSearchArgs and CWindowSizeArg set them
    // only then, blast_args.cpp:280-286,487-489) and `-no_greedy` (the algorithms); the
    // ungapped X-drop and the off-diagonal range are not recorded (`options.rs`).
    inputs.xdrop_gap = args.xdrop_gap;
    inputs.xdrop_gap_final = args.xdrop_gap_final;
    inputs.window_size = args.window_size;
    inputs.no_greedy = args.no_greedy;
    inputs.mb_index_name = args.megablast_index_name().map(str::to_string);
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
/// ```c++
///     SetDiagFilter(eDiagFilter_Post, "!(1305.31)");
///     CRef<CSeq_loc_Mapper> mapper = s_CreateSeqLocMapper(*query_id, seqloc,
///                                                         scope);
///     query_masks.Reset(mapper->Map(*query_masks));
///
///     const int kTopFlags = CSeq_loc::fStrand_Ignore|CSeq_loc::fMerge_All|CSeq_loc::fSort;
///     if (orig_query_mask.NotEmpty() && !orig_query_mask->IsNull()) {
///         CRef<CSeq_loc> tmp = orig_query_mask->Add(*query_masks,  kTopFlags, 0);
/// ```
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
/// ```c++
///         if (mqr) {
///             CRef<CSeqLocInfo> info_plus
///                 (new CSeqLocInfo(&* ival, CSeqLocInfo::eFramePlus1));
///             mqr->push_back(info_plus);
///
///             CRef<CSeqLocInfo> info_minus
///                 (new CSeqLocInfo(&* ival, CSeqLocInfo::eFrameMinus1));
///             mqr->push_back(info_minus);
///         }
/// ```
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
            let intervals = lcase_mask_intervals(
                record,
                &call.dust_masks[index],
                offset,
                call.lcase_masking,
                call.window_masker,
                call.both_strands,
                call.minus_strand,
            );
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{try_parse_from, Cli, Commands};
    use crate::report::archive::options::{self, algorithm_options, ParamValue};
    use crate::report::serial::real::double_to_string_posix;

    // The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote in the OF8 design's probes
    // (S/probe8/out, `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that
    // each archive adds: the probe's name, the options given besides the inputs, and the
    // list as `Name=variant:value`.
    const CASES: &[(&str, &[&str], &str)] = &[
    ("b_n_allN_db_fmt11", &[], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    // DI-6 (the frozen archives xdrop.*.11 of tests/outfmt0_manifest.tsv): the X-drops and the window
    // as given, -no_greedy's algorithms (eDynProgScoreOnly, eDynProgTbck = 0), the window of
    // dc-megablast replaced in its place, pyani ANIb's options.
    ("di6_nogreedy", &["-no_greedy", "-gapopen", "5", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:5 GapExtensionCost=integer:2 GapExtnAlgorithm=integer:0 GapTracebackAlgorithm=integer:0 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("di6_xdrops", &["-xdrop_ungap", "5", "-xdrop_gap", "12.5", "-xdrop_gap_final", "150", "-window_size", "40"], "EvalueThreshold=cutoff:e-value:10 GapXDropoff=real:12.5 GapXDropoffFinal=real:150 MaskAtHash=boolean:True UngappedMode=boolean:False WindowSize=integer:40 HitlistSize=integer:500"),
    ("di6_dc_win0", &["-task", "dc-megablast", "-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 MBTemplateType=integer:0 MBTemplateLength=integer:18 WordSize=integer:11 WindowSize=integer:0 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("di6_anib", &["-task", "blastn", "-xdrop_gap_final", "150", "-dust", "no", "-evalue", "1e-15", "-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:1e-15 GapXDropoffFinal=real:150 MaskAtHash=boolean:True DustFiltering=boolean:False UngappedMode=boolean:False HitlistSize=integer:1"),
    // DI-6b (the frozen archives minraw.*.11): the value is not recorded (the remote options'
    // eBlastOpt_CutoffScore is `if (0)`, api/blast_options_cxx.cpp:435-442).
    ("di6b_minraw", &["-min_raw_gapped_score", "15"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("e_n_badpair_fmt11", &["-reward", "2", "-penalty", "-7"], "EvalueThreshold=cutoff:e-value:10 MismatchPenalty=integer:-7 MatchReward=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("in_n_dust_minus_fmt11", &["-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("in_n_lcase_mask_fmt11", &["-lcase_masking"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    // SM-3 (the frozen archives sm3.*.11 of tests/outfmt0_manifest.tsv): the value as given.
    ("sm3_hard", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False UngappedMode=boolean:False HitlistSize=integer:500"),
    ("sm3_soft", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("in_n_mts1_fmt11", &["-task", "blastn", "-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:1"),
    ("in_n_qloc_fmt11", &["-query_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("in_n_qloc_full_fmt11", &["-query_loc", "1-480"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("m_n_dust_plus_fmt11", &["-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("m_n_lc_minus_fmt11", &["-lcase_masking", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g1", &["-word_size", "20", "-evalue", "1e-5", "-perc_identity", "90", "-max_hsps", "1", "-qcov_hsp_perc", "50"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:20 PercentIdentity=real:90 QueryCovHspPerc=real:50 MaxHspsPerSubject=integer:1 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g10", &["-gapopen", "2", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:2 GapExtensionCost=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g1rev", &["-qcov_hsp_perc", "50", "-max_hsps", "1", "-perc_identity", "90", "-evalue", "1e-5", "-word_size", "20"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:20 PercentIdentity=real:90 QueryCovHspPerc=real:50 MaxHspsPerSubject=integer:1 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g2", &["-penalty", "-3", "-culling_limit", "1", "-strand", "plus", "-dust", "no"], "EvalueThreshold=cutoff:e-value:10 MismatchPenalty=integer:-3 MaskAtHash=boolean:True DustFiltering=boolean:False UngappedMode=boolean:False Culling=integer:1 HitlistSize=integer:500"),
    ("opt_blastn_g3", &["-reward", "2", "-penalty", "-3", "-gapopen", "5", "-gapextend", "2", "-subject_besthit", "-lcase_masking"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:5 GapExtensionCost=integer:2 MismatchPenalty=integer:-3 MatchReward=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False UseSubjectBestHit=boolean:True HitlistSize=integer:500"),
    ("opt_blastn_g4", &["-max_target_seqs", "1", "-dust", "20 64 1", "-query_loc", "10-300", "-subject_loc", "10-300", "-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 MaskAtHash=boolean:True DustFilteringLevel=integer:20 DustFilteringWindow=integer:64 DustFilteringLinker=integer:1 UngappedMode=boolean:False HitlistSize=integer:1"),
    ("opt_blastn_g5", &["-num_descriptions", "5", "-num_alignments", "5", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:5"),
    ("opt_blastn_g6", &["-task", "dc-megablast", "-template_type", "coding", "-template_length", "16", "-word_size", "12"], "EvalueThreshold=cutoff:e-value:10 MBTemplateType=integer:0 MBTemplateLength=integer:16 WordSize=integer:12 WindowSize=integer:40 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g7", &["-task", "blastn-short"], "EvalueThreshold=cutoff:e-value:1000 MatchReward=integer:1 MismatchPenalty=integer:-3 WordSize=integer:7 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("short_max_hsps1", &["-task", "blastn-short", "-max_hsps", "1"], "EvalueThreshold=cutoff:e-value:1000 MatchReward=integer:1 MismatchPenalty=integer:-3 WordSize=integer:7 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True MaxHspsPerSubject=integer:1 UngappedMode=boolean:False HitlistSize=integer:500"),
    ("short_perc_identity90", &["-task", "blastn-short", "-perc_identity", "90"], "EvalueThreshold=cutoff:e-value:1000 MatchReward=integer:1 MismatchPenalty=integer:-3 WordSize=integer:7 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True PercentIdentity=real:90 UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g8", &["-task", "blastn", "-word_size", "7", "-dust", "yes"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:7 MaskAtHash=boolean:True DustFiltering=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g9", &["-task", "megablast", "-word_size", "28", "-evalue", "10", "-strand", "both", "-dust", "yes"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:28 MaskAtHash=boolean:True DustFiltering=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("optdb_blastn_g4", &["-max_target_seqs", "1", "-dust", "20 64 1", "-query_loc", "10-300", "-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 MaskAtHash=boolean:True DustFilteringLevel=integer:20 DustFilteringWindow=integer:64 DustFilteringLinker=integer:1 UngappedMode=boolean:False HitlistSize=integer:1"),
    ("t_n_cull0_fmt11", &["-culling_limit", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False Culling=integer:0 HitlistSize=integer:500"),
    ("t_n_dmb_opt_fmt11", &["-task", "dc-megablast", "-template_type", "optimal", "-template_length", "21"], "EvalueThreshold=cutoff:e-value:10 MBTemplateType=integer:1 MBTemplateLength=integer:21 WordSize=integer:11 WindowSize=integer:40 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("t_n_naln_fmt11", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:3"),
    ("t_n_ndesc_fmt11", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("t_n_pid0_fmt11", &["-perc_identity", "0"], "EvalueThreshold=cutoff:e-value:10 PercentIdentity=real:0 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("t_n_strandplus_lc_fmt11", &["-lcase_masking", "-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
];

    fn posix(value: f64) -> String {
        let mut digits = Vec::new();
        double_to_string_posix(value, 15, &mut digits);
        String::from_utf8(digits).unwrap()
    }

    fn render(params: &[crate::report::archive::options::Param]) -> String {
        params
            .iter()
            .map(|param| {
                let value = match &param.value {
                    ParamValue::Evalue(value) => format!("cutoff:e-value:{}", posix(*value)),
                    ParamValue::Integer(value) => format!("integer:{value}"),
                    ParamValue::BigInteger(value) => format!("big-integer:{value}"),
                    ParamValue::Boolean(true) => "boolean:True".to_string(),
                    ParamValue::Boolean(false) => "boolean:False".to_string(),
                    ParamValue::Real(value) => format!("real:{}", posix(*value)),
                    ParamValue::Str(value) => format!("string:{value}"),
                    ParamValue::IntegerList(values) => format!("integer-list:{values:?}"),
                };
                format!("{}={value}", param.name())
            })
            .collect::<Vec<_>>()
            .join(" ")
    }

    fn parse(extra: &[&str]) -> Option<BlastnArgs> {
        let mut argv = vec!["losat", "blastn", "-query", "q.fa", "-subject", "s.fa"];
        argv.extend_from_slice(extra);
        argv.extend_from_slice(&["-outfmt", "11"]);
        let cli: Cli = try_parse_from(argv).ok()?;
        match cli.command {
            Commands::Blastn(args) => Some(args),
            _ => None,
        }
    }

    #[test]
    fn the_option_lists_equal_those_of_ncbis_probes() {
        let mut compared = 0;
        for &(name, extra, expected) in CASES {
            let args = parse(extra).unwrap_or_else(|| panic!("{name} {extra:?} does not parse"));
            let dust = args.dust_filtering.clone();
            let inputs = option_inputs(&args, dust.as_deref());
            assert_eq!(
                render(&algorithm_options(&inputs)),
                expected,
                "{name} {extra:?}"
            );
            compared += 1;
        }
        assert!(compared >= 27, "{compared}");
    }

    fn dust(ranges: &[(usize, usize)]) -> Vec<MaskedInterval> {
        ranges
            .iter()
            .map(|&(start, end)| MaskedInterval::new(start, end))
            .collect()
    }

    #[test]
    fn lcase_mask_is_the_query_vector_mask() {
        use crate::blastinput::fasta_reader::RecordSource;
        // Lower case at 0..20 of the record, the range starts at 29: the lower case outside
        // the range is kept, in record coordinates, once per strand.
        let mut whole = vec![b'a'; 20];
        whole.extend(std::iter::repeat_n(b'C', 280));
        let mut record = FastaRecord::new("Query_1", b"q", &whole[29..]);
        record.source = Some(Box::new(RecordSource {
            defline: b">q".to_vec(),
            uncut: Some(whole),
            ..RecordSource::default()
        }));
        assert_eq!(
            lcase_mask_intervals(&record, &[], 29, true, false, true, false),
            vec![(0, 19), (0, 19)]
        );
        // Without -lcase_masking there is no user mask.
        assert!(lcase_mask_intervals(&record, &[], 29, false, false, true, false).is_empty());
        // DUST intervals of the range are mapped by the range's start, merged with the lower
        // case when they overlap or abut, and sorted.
        assert_eq!(
            lcase_mask_intervals(
                &record,
                &dust(&[(0, 5), (30, 40)]),
                20,
                true,
                false,
                true,
                false
            ),
            vec![(0, 24), (0, 24), (50, 59), (50, 59)]
        );
        // A DUST-only mask: ascending once on a single strand, descending with -strand minus,
        // twice with both strands.
        let plain = FastaRecord::new("Query_2", b"q", &[b'A'; 100]);
        let intervals = dust(&[(5, 10), (20, 30)]);
        assert_eq!(
            lcase_mask_intervals(&plain, &intervals, 0, false, false, false, false),
            vec![(5, 9), (20, 29)]
        );
        assert_eq!(
            lcase_mask_intervals(&plain, &intervals, 0, false, false, false, true),
            vec![(20, 29), (5, 9)]
        );
        assert_eq!(
            lcase_mask_intervals(&plain, &intervals, 0, false, false, true, false),
            vec![(5, 9), (5, 9), (20, 29), (20, 29)]
        );
        // A query with lower case is ascending even with -strand minus.
        let mut lower = vec![b'A'; 100];
        lower[50] = b'c';
        let lower = FastaRecord::new("Query_3", b"q", &lower);
        assert_eq!(
            lcase_mask_intervals(&lower, &intervals, 0, true, false, false, true),
            vec![(5, 9), (5, 9), (20, 29), (20, 29), (50, 50), (50, 50)]
        );
        // WindowMasker's masks replace the others: each interval twice (the plus and minus
        // frames), in the masker's order, with every strand option (NCBI 2.17.0's archives of
        // `-window_masker_db`, line NP M2).
        for (both, minus) in [(true, false), (false, false), (false, true)] {
            assert_eq!(
                lcase_mask_intervals(&lower, &intervals, 3, true, true, both, minus),
                vec![(8, 12), (8, 12), (23, 32), (23, 32)]
            );
        }
    }

    #[test]
    fn the_service_follows_the_task() {
        for (task, service) in [
            ("megablast", "megablast"),
            ("dc-megablast", "megablast"),
            ("blastn", "plain"),
            ("blastn-short", "plain"),
        ] {
            assert_eq!(options::service(Program::Blastn, task), service);
        }
    }

    #[test]
    fn error_lines_are_ncbis_error_texts_only() {
        let native = |message: &str| {
            anyhow::Error::from(crate::cli::NativeError {
                exit: 1,
                message: message.to_string().into(),
            })
        };
        let lines = error_lines(
            &native(
                "BLAST query/options error: bad  option\nPlease refer to the LOSAT user manual.\n",
            ),
            false,
        )
        .unwrap();
        assert_eq!(
            lines,
            vec![
                (ERROR, b"BLAST query/options error: bad option".to_vec()),
                (ERROR, b"Please refer to the LOSAT user manual.".to_vec()),
            ]
        );
        assert!(error_lines(&native("Command line argument error: x\n"), false).is_none());
        let subject = error_lines(
            &native("Command line argument error: Argument \"subject\". File is not accessible:  `a.fa'\n"),
            true,
        )
        .unwrap();
        assert_eq!(
            subject,
            vec![(
                ERROR,
                b"Command line argument error: Argument \"subject\". File is not accessible: `a.fa'".to_vec()
            )]
        );
        assert!(error_lines(&anyhow::anyhow!("unsupported by LOSAT"), false).is_none());
        // `CATCH_ALL`'s report of a `std::exception` (exit status 255, blast_app_util.hpp:260-263):
        // one error message (line RG D3-1, `Error: vector::reserve` of the `-in_pssm` reader).
        let std_exception =
            anyhow::Error::from(crate::cli::NativeError::std_exception("vector::reserve"));
        assert_eq!(
            error_lines(&std_exception, false).unwrap(),
            vec![(ERROR, b"Error: vector::reserve".to_vec())]
        );
        // A text with `Error: ` and another exit status is not one of `CATCH_ALL`'s reports.
        assert!(error_lines(&native("Error: something\n"), false).is_none());
    }
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
