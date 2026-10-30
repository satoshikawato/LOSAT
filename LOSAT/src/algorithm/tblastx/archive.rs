//! TBLASTX's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:198
//! ```c++
//!                 if (fmt_args->ArchiveFormatRequested(args)) {
//!                     formatter.WriteArchive(*queries, *opts_hndl, *results, 0, m_Bah.GetMessages());
//!                     m_Bah.ResetMessages();
//!                 } else {
//! ```
//! One archive per query batch (the reader's batch, `BATCH_SIZE` letters, 10002 by default),
//! written after the batch's search. LOSAT's TBLASTX searches the batches first and writes every
//! format at the end, from the hits of all of them (`write_tblastx_outputs`), so the archives are
//! written there, one per batch that was read and searched, with that batch's searched flag.
//! `report::archive` writes the document and this module gathers what TBLASTX holds for its
//! other formats: the Seq-aligns of outfmt 8, the SEG masks of the queries (in nucleotide
//! coordinates, one mask per frame as BLASTX), and the Karlin-Altschul block and search space of
//! the batch's first query.
//!
//! The batch-wide values: a batch is searched on its own (`search_query_batch` takes the
//! batch's queries), so the values NCBI computes over its batch (the context lengths and
//! Karlin blocks, the average query length and the smallest lambda of the linking cutoffs, the
//! cutoffs) are the batch's, as for outfmt 8. The archive reads the batch's first query only:
//! its first valid context's Karlin block, search space and length adjustment
//! (`TblastxQueryStats`), and the effective search space of context 0 (the first query's frame
//! +1: 0 when that context is not valid), which is the `EffectiveSearchSpace` option of the
//! batch (`SplitQuery_SetEffectiveSearchSpace`, `BlastOptions` memento: the first context).

use std::borrow::Cow;
use std::io::Write;
use std::ops::Range;

use anyhow::Result;

use super::blast_engine::{TblastxHsp, TblastxQueryStats};
use super::report::SubjectLetters;
use super::seqalign::{query_records_of, SeqalignSubjects};
use super::subjects::SubjectSetNumbers;
use super::TblastxArgs;
use crate::algorithm::tblastn::archive::lowercase_runs;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::blastinput::value_parsers::SegSpec;
use crate::report::archive::bioseq::BioseqInput;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{Filter, OptionInputs};
use crate::report::archive::{
    defline_of, ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive,
    KaBlock,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;
use crate::report::serial::seqalign::{write_seq_aligns, Segments};

/// The recorded options of a TBLASTX run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2). TBLASTX has one task, and the options it does not
/// accept (`-strand`, ...) are never recorded.
pub(crate) fn option_inputs(args: &TblastxArgs) -> OptionInputs {
    let mut inputs = OptionInputs::new(Program::Tblastx, "tblastx");
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    // `-evalue` replaces the handle's value in place: the same list whether it was given.
    inputs.evalue = Some(args.evalue);
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps.map(|count| count as i32);
    // `-matrix` is recorded as typed when given; `-threshold` is read as a double and recorded as
    // an `int` (`static_cast<int>`). Without `-threshold` and `-window_size` the list holds the
    // matrix's suggestions (`report::archive::options`).
    inputs.matrix = args.matrix.clone();
    inputs.threshold = args.threshold.map(|threshold| threshold as i32);
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
    // ```c++
    //     if (args[kArgLookupTableMaskingOnly]) {
    //         opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
    //     }
    // ```
    // `MaskAtHash` is recorded with the value given (the class always sets it: false when absent).
    inputs.soft_masking = args.soft_masking;
    inputs.window_size = args.window_size;
    inputs.seg = args.seg_spec().ok().map(|spec| match spec {
        SegSpec::No => Filter::No,
        SegSpec::Yes => Filter::Yes,
        SegSpec::WindowLocutHicut {
            window,
            locut,
            hicut,
        } => Filter::Seg {
            window,
            locut,
            hicut,
        },
    });
    inputs.culling_limit = args.culling_limit;
    inputs.query_gencode = i32::from(args.query_gencode);
    inputs.db_gencode = i32::from(args.db_gencode);
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(500, |size| size as i32);
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with TBLASTX's recorded options).
pub(crate) fn handle_for(
    args: &TblastxArgs,
    archive_requested: bool,
) -> Result<Option<ArchiveHandle>> {
    if !archive_requested {
        return Ok(None);
    }
    // NCBI keeps the name without the white space at its ends (`NStr::TruncateSpaces`,
    // blast_args.cpp:2472-2476), and the request's `database` is that name.
    ArchiveHandle::for_program(&option_inputs(args), args.database_name(), "TBLASTX").map(Some)
}

/// `report::archive::check_titles` for TBLASTX's records.
pub(crate) fn check_titles(records: &[FastaRecord], role: &str) -> Result<()> {
    crate::report::archive::check_titles(records, role, "TBLASTX")
}

pub(crate) use crate::report::archive::error_lines;

fn frame_of(frame: i8) -> Frame {
    match frame {
        1 => Frame::Plus1,
        2 => Frame::Plus2,
        3 => Frame::Plus3,
        -1 => Frame::Minus1,
        -2 => Frame::Minus2,
        _ => Frame::Minus3,
    }
}

/// NCBI's order of the frames of a query's masks (the contexts: plus1, plus2, plus3, minus1,
/// minus2, minus3).
fn frame_rank(frame: i8) -> u8 {
    if frame > 0 {
        frame as u8
    } else {
        3 + frame.unsigned_abs()
    }
}

/// The results masks of a searched query: one mask per frame with intervals, in NCBI's frame
/// order, the intervals in nucleotide coordinates of the searched letters moved into the record
/// by the start of the `-query_loc` interval; a mask that covers the whole query is dropped
/// (`range != kTarget`, blast_aux.cpp:936-951), as BLASTX's.
fn query_results_masks(
    query: &FastaRecord,
    stats: &TblastxQueryStats,
    offset: usize,
) -> Vec<MaskEntry> {
    let length = query.seq().len() as i32;
    let mut frames = super::report::query_dna_masks(stats, query.seq().len());
    frames.sort_by_key(|(frame, _)| frame_rank(*frame));
    let offset = offset as i32;
    frames
        .into_iter()
        .filter_map(|(frame, intervals)| {
            let intervals: Vec<(i32, i32)> = intervals
                .into_iter()
                .filter(|&(from, to)| from <= to && !(from == 0 && to == length - 1))
                .map(|(from, to)| (from + offset, to + offset))
                .collect();
            (!intervals.is_empty()).then(|| MaskEntry {
                id: query.local_id.clone(),
                frame: frame_of(frame),
                intervals,
            })
        })
        .collect()
}

/// The report data of a run, as `write_tblastx_outputs` holds it.
pub(super) struct ArchiveData<'a> {
    /// The searched queries of every batch written, in order (each cut to its `-query_loc`).
    pub queries: &'a [FastaRecord],
    pub query_stats: &'a [TblastxQueryStats],
    /// Whether each query's batch was searched.
    pub searched: &'a [bool],
    /// `-subject`: the records (every archive writes all of them whole); empty for `-db`.
    pub subject_records: &'a [FastaRecord],
    /// The final hits in report order and their engine counts.
    pub hits: &'a [TblastxHsp],
    pub engine: &'a [(u32, u32)],
    pub subjects: &'a SeqalignSubjects<'a>,
    pub subject_letters: SubjectLetters<'a>,
    pub query_placements: &'a Placements,
    pub numbers: &'a SubjectSetNumbers,
    /// The queries of each batch written, in order.
    pub spans: &'a [Range<usize>],
    /// `-lcase_masking`.
    pub lcase_masking: bool,
}

/// The Seq-aligns of one batch: the elements of outfmt 8, one query at a time.
struct TblastxAligns<'a> {
    run: RunRecord<'a>,
    data: &'a ArchiveData<'a>,
    span: Range<usize>,
}

impl<'a> AlignSource for TblastxAligns<'a> {
    fn write<E: crate::report::serial::encoder::Encoder>(
        &mut self,
        out: &mut ArchiveOut<'_, E>,
    ) -> Result<()> {
        let data = self.data;
        let mut segments = Segments::default();
        let mut subjects: Vec<SubjectRecord<'a>> = Vec::new();
        let mut hsps: Vec<HspRecord<'a>> = Vec::new();
        for q_idx in self.span.clone() {
            // The final list is in query order: this query's run of hits.
            let from = data
                .hits
                .partition_point(|hsp| (hsp.hit.q_idx as usize) < q_idx);
            let to =
                from + data.hits[from..].partition_point(|hsp| (hsp.hit.q_idx as usize) == q_idx);
            if from == to {
                continue;
            }
            query_records_of(
                &data.hits[from..to],
                &data.engine[from..to],
                data.subjects,
                data.subject_letters,
                &mut subjects,
                &mut hsps,
            )?;
            let query = &data.queries[q_idx];
            let query = QueryRecord {
                id: SeqId::Local(&query.local_id),
                searched_length: i32::try_from(query.seq().len())?,
                offset: i32::try_from(data.query_placements.offset(q_idx))?,
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

/// Writes the archive of each batch of `data.spans` to `writer`; `warnings` writes the warnings
/// of a batch's queries (to the diagnostics, which saves them as the batch's messages) before
/// its archive.
pub(super) fn write_tblastx_archives(
    handle: &ArchiveHandle,
    data: &ArchiveData<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    writer: &mut dyn Write,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Tblastx,
        db_title: match data.subjects {
            SeqalignSubjects::Records { .. } => None,
            SeqalignSubjects::Database { title, .. } => Some(title),
        },
        gapped: false,
    };
    // The Bioseq inputs borrow their deflines: keep them alive for the whole call.
    let subject_deflines: Vec<Cow<'_, [u8]>> =
        data.subject_records.iter().map(defline_of).collect();
    let subjects_input: Vec<BioseqInput<'_>> = data
        .subject_records
        .iter()
        .zip(&subject_deflines)
        .map(|(record, defline)| BioseqInput {
            id: &record.local_id,
            defline,
            title: &record.title,
            letters: record.archive_letters(),
            protein: false,
        })
        .collect();
    for span in data.spans.iter().cloned() {
        if span.is_empty() {
            continue;
        }
        if let Some(warnings) = warnings.as_deref_mut() {
            for index in span.clone() {
                warnings.before_query(index, &mut *writer)?;
            }
        }
        let records = &data.queries[span.clone()];
        let deflines: Vec<Cow<'_, [u8]>> = records.iter().map(defline_of).collect();
        let queries: Vec<BioseqInput<'_>> = records
            .iter()
            .zip(&deflines)
            .map(|(record, defline)| BioseqInput {
                id: &record.local_id,
                defline,
                title: &record.title,
                letters: record.archive_letters(),
                protein: false,
            })
            .collect();
        let first = span.start;
        let searched = data.searched[first];

        // `RequiredStart` / `RequiredEnd`: a `-query_loc` that is not the first query's whole
        // length.
        let placements = data.query_placements;
        let required = if placements.ranged() {
            let placement = &placements.0[first];
            let cut = data.queries[first].seq().len();
            (cut != placement.length)
                .then_some((placement.offset as i32, (placement.offset + cut) as i32 - 1))
        } else {
            None
        };

        // The results' masks are those of the search, so only a searched batch has them.
        let mut results_masks: Vec<MaskEntry> = Vec::new();
        if searched {
            for (index, record) in (span.start..).zip(records) {
                results_masks.extend(query_results_masks(
                    record,
                    &data.query_stats[index],
                    placements.offset(index),
                ));
            }
        }

        // The batch's first query: its first valid context gives the block, the search space
        // and the length adjustment (a first query without one gives both blocks the `-1`
        // triple); the `EffectiveSearchSpace` option is the space of context 0 (frame +1).
        let stats = &data.query_stats[first];
        let (ka_blocks, search_space, length_adjustment) = match &stats.karlin {
            Some(karlin) => {
                let adjustment = data
                    .numbers
                    .context_lengths(stats.aa_len as i64, karlin)
                    .map_or(0, |lengths| lengths.length_adjustment);
                (
                    [ka_block(karlin, false), KaBlock::missing(true)],
                    stats.eff_searchsp,
                    adjustment,
                )
            }
            None => ([KaBlock::missing(false), KaBlock::missing(true)], 0, 0),
        };
        // `-lcase_masking`: the `LCaseMask` parameters of the request are the lower-case runs of
        // each whole query record, each one for both strands (as BLASTX's `lcase_masks`).
        let mut lcase_masks: Vec<MaskEntry> = Vec::new();
        if data.lcase_masking {
            for record in records {
                for run in lowercase_runs(record.archive_letters()) {
                    for frame in [Frame::Plus1, Frame::Minus1] {
                        lcase_masks.push(MaskEntry {
                            id: record.local_id.clone(),
                            frame,
                            intervals: vec![run],
                        });
                    }
                }
            }
        }
        let batch = BatchArchive {
            queries: &queries,
            subjects: &subjects_input,
            required,
            lcase_masks: &lcase_masks,
            masks: &results_masks,
            ka_blocks,
            search_space,
            length_adjustment,
            searched_space: searched.then_some(stats.context0_search_space),
        };
        let mut aligns = TblastxAligns {
            run: RunRecord {
                program: run.program,
                db_title: run.db_title,
                gapped: run.gapped,
            },
            data,
            span,
        };
        write_archive(handle, &batch, &mut aligns, writer)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{try_parse_from, Cli, Commands};
    use crate::report::archive::options::{algorithm_options, ParamValue};
    use crate::report::serial::real::double_to_string_posix;

    // The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote for tblastx (S/of8/b3/out, the OF8b3
    // probes; `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that each
    // archive adds: the name, the options given besides the inputs, and the list as
    // `Name=variant:value`. LOSAT's TBLASTX refuses some of these options at the command line
    // (`-sum_stats`, `-strand`, `-subject_besthit`, `-xdrop_ungap`,
    // `-parse_deflines`) or in its limits (`-word_size` other than 3): the list that NCBI would
    // record is what the archive's model knows about them, pinned here.
    const CASES: &[(&str, &[&str], &str)] = &[
        ("def_sub", &[], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("seg_no_sub", &["-seg", "no"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("seg_yes_sub", &["-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:True WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("seg_par_sub", &["-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("seg_def_sub", &[], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("batch100_sub", &[], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("dgc11_sub", &["-db_gencode", "11"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:11 HitlistSize=integer:500"),
        ("dgc2_db", &["-db_gencode", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:2 HitlistSize=integer:500"),
        ("dgc4_db", &["-db_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:4 HitlistSize=integer:500"),
        ("gc44_sub", &["-query_gencode", "4", "-db_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:4 DbGeneticCode=integer:4 HitlistSize=integer:500"),
        ("mts_warn_sub", &["-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:2"),
        ("opt_besthit_sub", &["-subject_besthit"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 UseSubjectBestHit=boolean:True WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_cull0_sub", &["-culling_limit", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:0 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_cull2_db", &["-culling_limit", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:2 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_cull5_sub", &["-culling_limit", "5"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:5 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_cull_mts_sub", &["-culling_limit", "3", "-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:3 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:2"),
        ("opt_dbsize_db", &["-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_ev10_sub", &["-evalue", "10"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_ev_sub", &["-evalue", "0.01"], "EvalueThreshold=cutoff:e-value:0.01 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_evbig_sub", &["-evalue", "1000"], "EvalueThreshold=cutoff:e-value:1000 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_hsps_sub", &["-max_hsps", "1"], "EvalueThreshold=cutoff:e-value:10 MaxHspsPerSubject=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_m45_sub", &["-matrix", "BLOSUM45"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM45 WordThreshold=integer:16 WindowSize=integer:60 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_m62_sub", &["-matrix", "BLOSUM62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM62 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_m62lc_sub", &["-matrix", "blosum62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:blosum62 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_m80_sub", &["-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WordThreshold=integer:14 WindowSize=integer:25 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_many_db", &["-evalue", "1e-5", "-threshold", "12", "-window_size", "30", "-max_hsps", "2", "-qcov_hsp_perc", "30", "-query_gencode", "4", "-db_gencode", "4", "-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:1e-05 QueryCovHspPerc=real:30 MaxHspsPerSubject=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:12 WindowSize=integer:30 QueryGeneticCode=integer:4 DbGeneticCode=integer:4 HitlistSize=integer:500"),
        ("opt_many_sub", &["-evalue", "1e-5", "-word_size", "3", "-threshold", "12", "-window_size", "30", "-max_hsps", "2", "-qcov_hsp_perc", "30", "-query_gencode", "4", "-db_gencode", "4", "-seg", "12 2.2 2.5", "-strand", "plus", "-culling_limit", "4"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:3 QueryCovHspPerc=real:30 MaxHspsPerSubject=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 Culling=integer:4 WindowSize=integer:30 QueryGeneticCode=integer:4 DbGeneticCode=integer:4 HitlistSize=integer:500"),
        ("opt_mts1_db", &["-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:1"),
        ("opt_mts4_sub", &["-max_target_seqs", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:4"),
        ("opt_naln_sub", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:3"),
        ("opt_ndesc_sub", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_ndescaln_sub", &["-num_descriptions", "3", "-num_alignments", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:4"),
        ("opt_pam30_sub", &["-matrix", "PAM30"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM30 WordThreshold=integer:18 WindowSize=integer:15 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_pam70_sub", &["-matrix", "PAM70"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM70 WordThreshold=integer:16 WindowSize=integer:20 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_parse_sub", &["-parse_deflines"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_qcov_sub", &["-qcov_hsp_perc", "50"], "EvalueThreshold=cutoff:e-value:10 QueryCovHspPerc=real:50 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_soft_sub", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:True SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_softf_sub", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_sumstats_sub", &["-sum_stats", "true"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:True LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_sumstatsf_sub", &["-sum_stats", "false"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_thr13_sub", &["-threshold", "13"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_thr_frac_sub", &["-threshold", "11.7"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:11 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_thr_sub", &["-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:15 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_win0_sub", &["-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:0 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_win40_sub", &["-window_size", "40"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_win_sub", &["-window_size", "30"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:30 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_ws2_sub", &["-word_size", "2"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_ws3_db", &["-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_ws4_sub", &["-word_size", "4"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:4 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("opt_xdropu_sub", &["-xdrop_ungap", "20"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("qgc11_sub", &["-query_gencode", "11"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:11 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("qgc2_sub", &["-query_gencode", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:2 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("qgc4_db", &["-query_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:4 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("qloc_db", &["-query_loc", "10-200"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("qloc_full_sub", &["-query_loc", "1-366"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("qloc_minus_sub", &["-query_loc", "10-200", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("qsloc_sub", &["-query_loc", "10-200", "-subject_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("sloc_sub", &["-subject_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("strand_both_sub", &["-strand", "both"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("strand_minus_db", &["-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
        ("strand_plus_db", &["-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
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
                };
                format!("{}={value}", param.name)
            })
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// The recorded inputs of a command line. The options that LOSAT's TBLASTX does not accept
    /// (they stop at the command line) are applied to the inputs by hand.
    fn inputs_of(extra: &[&str]) -> Option<OptionInputs> {
        let mut parsed: Vec<&str> = Vec::new();
        let mut by_hand: Vec<(&str, &str)> = Vec::new();
        let mut index = 0;
        while index < extra.len() {
            match extra[index] {
                "-sum_stats" => {
                    by_hand.push((extra[index], extra[index + 1]));
                    index += 2;
                }
                "-subject_besthit" => {
                    by_hand.push((extra[index], ""));
                    index += 1;
                }
                // Not recorded (the ungapped X-drop is never set by the application, the strand
                // is not an option of the request).
                "-xdrop_ungap" | "-strand" => index += 2,
                "-parse_deflines" => index += 1,
                other => {
                    parsed.push(other);
                    index += 1;
                }
            }
        }
        let mut argv = vec!["losat", "tblastx", "-query", "q.fa", "-subject", "s.fa"];
        argv.extend_from_slice(&parsed);
        argv.extend_from_slice(&["-outfmt", "11"]);
        let cli: Cli = try_parse_from(argv).ok()?;
        let Commands::Tblastx(args) = cli.command else {
            return None;
        };
        let mut inputs = option_inputs(&args);
        for (name, value) in by_hand {
            match name {
                "-sum_stats" => inputs.sum_stats = Some(value == "true"),
                _ => inputs.subject_besthit = true,
            }
        }
        Some(inputs)
    }

    #[test]
    fn the_option_lists_equal_those_of_ncbis_probes() {
        assert!(CASES.len() >= 55);
        for &(name, extra, expected) in CASES {
            let inputs =
                inputs_of(extra).unwrap_or_else(|| panic!("{name} {extra:?} does not parse"));
            assert_eq!(
                render(&algorithm_options(&inputs)),
                expected,
                "{name} {extra:?}"
            );
        }
    }

    #[test]
    fn the_service_is_plain() {
        assert_eq!(
            crate::report::archive::options::service(Program::Tblastx, "tblastx"),
            "plain"
        );
    }

    #[test]
    fn results_masks_are_per_frame_in_ncbis_order_without_a_whole_query_mask() {
        // A 30-nt query: frame lengths 10, 9, 9 (BLAST_GetTranslatedProteinLength), so the
        // DNA length is 30; frames come in the order of the contexts, not of the stats.
        let record = FastaRecord::new("Query_1", b"q", &[b'A'; 30]);
        let stats = TblastxQueryStats {
            seg_masks: vec![(-1, vec![(2, 5)]), (1, vec![(2, 5)]), (2, vec![(0, 9)])],
            ..TblastxQueryStats::default()
        };
        let masks = query_results_masks(&record, &stats, 0);
        let frames: Vec<Frame> = masks.iter().map(|mask| mask.frame).collect();
        assert_eq!(frames, [Frame::Plus1, Frame::Plus2, Frame::Minus1]);
        // Residues 2..=4 of frame +1: NCBI's end is 3 * 4 (blast_filter.c:925-931), the formula
        // `query_masks_convert_with_ncbis_frame_formulas` pins.
        assert_eq!(masks[0].intervals, [(6, 12)]);
        // A mask that spans the whole query is dropped (blast_aux.cpp:936-951), and a
        // `-query_loc` start moves the rest into the record.
        let moved = query_results_masks(&record, &stats, 100);
        assert_eq!(moved[0].intervals, [(106, 112)]);
    }

    #[test]
    fn frames_rank_as_ncbis_contexts() {
        let ranks: Vec<u8> = [1i8, 2, 3, -1, -2, -3]
            .into_iter()
            .map(frame_rank)
            .collect();
        assert_eq!(ranks, [1, 2, 3, 4, 5, 6]);
    }
}
