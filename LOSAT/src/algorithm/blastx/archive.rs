//! BLASTX's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:284-286
//! ```c++
//!             if (fmt_args->ArchiveFormatRequested(args)) {
//!                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//!                 m_Bah.ResetMessages();
//!             } else {
//! ```
//! One archive per query batch (the reader's batch, `BATCH_SIZE` letters), written after the
//! batch's search; `report::archive` writes the document, this module gathers what BLASTX holds
//! for its other formats: the hit lists (the Seq-aligns of outfmt 8), the masks of the six
//! translated frames (`QueryContext::dna_masks`, as outfmt 15 and 16 print them), and the
//! Karlin-Altschul blocks and search space of the batch's first query. A query has six
//! contexts (frames +1, +2, +3, -1, -2, -3); the batch's `EffectiveSearchSpace` option is the
//! search space of the first query's context 0 (frame +1; 0 when that context is not searched:
//! `-strand minus`, a query shorter than a codon, a frame the Karlin-Altschul statistics reject),
//! while `search-stats` and the blocks are those of its first valid context.

use std::borrow::Cow;
use std::io::Write;

use anyhow::Result;

use super::args::{BlastxArgs, ResolvedOptions};
use super::input::FastaRecord;
use super::results::HspList;
use super::runtime::BatchResults;
use super::seqalign::{query_records, SeqalignSubjects};
use crate::blastinput::app::Program;
use crate::report::archive::bioseq::BioseqInput;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{Filter, OptionInputs};
use crate::report::archive::{
    ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive, KaBlock,
};
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;
use crate::report::serial::seqalign::{write_seq_aligns, Segments};
use crate::stats::tables::KarlinParams;

/// The recorded options of a BLASTX run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2).
pub(crate) fn option_inputs(args: &BlastxArgs) -> OptionInputs {
    let mut inputs = OptionInputs::new(Program::Blastx, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    inputs.evalue = Some(args.evalue);
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.matrix = args.matrix.clone();
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps;
    inputs.sum_stats = args.sum_stats;
    // `-threshold` is read as a double and recorded as an `int` (`static_cast<int>`).
    inputs.threshold = args.threshold.map(|threshold| threshold as i32);
    inputs.window_size = args.window_size;
    inputs.ungapped = args.ungapped;
    inputs.soft_masking = Some(args.soft_masking);
    inputs.seg = Some(if !args.seg.enabled {
        Filter::No
    } else if args.seg.yes {
        Filter::Yes
    } else {
        Filter::Seg {
            window: args.seg.window,
            locut: args.seg.locut,
            hicut: args.seg.hicut,
        }
    });
    inputs.culling_limit = args.culling_limit;
    inputs.subject_besthit = args.subject_besthit;
    inputs.max_intron_length = args.max_intron_length;
    inputs.query_gencode = i32::from(args.query_gencode);
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments.map(|size| size as i32))
        .map_or(500, |size| size);
    inputs.comp_based_stats = i32::from(args.comp_based_stats);
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with BLASTX's recorded options).
pub(crate) fn handle_for(
    args: &BlastxArgs,
    archive_requested: bool,
) -> Result<Option<ArchiveHandle>> {
    if !archive_requested {
        return Ok(None);
    }
    // NCBI keeps the name without the white space at its ends (`NStr::TruncateSpaces`,
    // blast_args.cpp:2472-2476), and the request's `database` is that name.
    let database = args.db.as_deref().map(crate::seqdb::db_argument_name);
    ArchiveHandle::for_program(&option_inputs(args), database, "BLASTX").map(Some)
}

/// The refusal of a title that the archive cannot write as NCBI does
/// (`report::archive::ArchiveEncoding::writes_text`: in text and binary a byte outside
/// `0x20..=0x7E`, in XML a control byte) for BLASTX's records, before the first byte.
pub(crate) fn check_titles(records: &[FastaRecord], role: &str) -> Result<()> {
    let encoding = crate::report::archive::ArchiveEncoding::from_env();
    if let Some(record) = records
        .iter()
        .find(|record| !encoding.writes_text(record.title.as_bytes()))
    {
        anyhow::bail!(
            "the title of the {role} {} has a character the archive cannot write as NCBI does (outside the visible ASCII range in text and binary, a control character in XML), which is not supported by LOSAT's BLASTX with -outfmt 11",
            record.internal_id
        );
    }
    Ok(())
}

pub(crate) use crate::report::archive::error_lines;

/// A title the archive cannot write, found by `prescan_query_titles` (`check_titles`); a typed
/// error, so that it is told apart from the reader's own errors, which end the scan.
#[derive(Debug)]
struct TitleRefusal(String);

impl std::fmt::Display for TitleRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TitleRefusal {}

/// Whether the query is a regular file, decided on the handle already open (`fstat`), before any
/// second `open` by path: `prescan_query_titles` reads the file a second time, which would take
/// bytes from a pipe (`<(...)`, `/dev/stdin`) or wait for a writer of a named FIFO.
pub(crate) fn query_is_regular(file: &std::fs::File) -> bool {
    file.metadata().is_ok_and(|metadata| metadata.is_file())
        && crate::blastinput::input_files::has_stream_position(file)
}

/// The refusal of a query for which `query_is_regular` is false.
pub(crate) fn refuse_irregular_query() -> Result<()> {
    anyhow::bail!("unsupported BLASTX query input for -outfmt 11: use a regular file")
}

/// Reads the query file once before the search, batch by batch as the search reads it
/// (`batch_size` letters, the run's), and refuses the first record whose title the archive
/// cannot write (`check_titles`), so that no refusal can follow a written byte. A reader error
/// ends the scan: the run reads the file again and fails there as it always did, after the
/// archives of the batches before it (NCBI writes them), so the titles of the batches the reader
/// completed before the error are checked; the records of the batch the error interrupts are
/// never written and need no check. The file must be a regular file (`query_is_regular`).
pub(crate) fn prescan_query_titles(
    path: &std::path::Path,
    lowercase: bool,
    batch_size: usize,
) -> Result<()> {
    let file = std::fs::File::open(path)?;
    let mut stream = super::input::FastaStream::from_file(file);
    let scanned = super::input::parse_fasta_stream_with_warnings(
        &mut stream,
        false,
        lowercase,
        batch_size,
        &mut |batch| {
            check_titles(batch, "query").map_err(|refusal| TitleRefusal(refusal.to_string()).into())
        },
        &mut |_| Ok(()),
    );
    match scanned {
        Err(error) if error.is::<TitleRefusal>() => Err(error),
        _ => Ok(()),
    }
}

/// What the archive of one query batch needs.
pub(crate) struct ArchiveCall<'a> {
    pub handle: &'a ArchiveHandle,
    pub options: &'a ResolvedOptions,
    /// The batch's queries; `BatchResults::query_ordinal` indexes them.
    pub queries: &'a [FastaRecord],
    /// The subjects by the oids of the hit lists.
    pub subjects: &'a [FastaRecord],
    pub names: SeqalignSubjects<'a>,
    /// `-subject`: every subject record (the Bioseqs of the request); `None` for `-db`.
    pub subject_records: Option<&'a [FastaRecord]>,
    /// The results of the batch: one entry, or several when the search split a batch above
    /// 10002 letters (`BATCH_SIZE`).
    pub batches: &'a [BatchResults],
}

/// The values of a batch's first query that `search-stats`, `ka-blocks` and the
/// `EffectiveSearchSpace` option take.
struct FirstContext {
    /// The ungapped block of the first valid context and the gapped one (absent for an
    /// ungapped search); `None` when the first query has no valid context (both blocks are
    /// then missing).
    karlin: Option<(KarlinParams, Option<KarlinParams>)>,
    /// The effective search space and length adjustment of that context (0 without one).
    search_space: i64,
    length_adjustment: i64,
    /// The search space of the first query's context 0 (the option); 0 when it is not valid.
    context0_space: i64,
}

fn first_context(batch: &BatchResults, gapped: Option<KarlinParams>) -> FirstContext {
    let parameters = &batch.parameters[..6];
    let valid = parameters.iter().find(|parameters| parameters.valid);
    FirstContext {
        karlin: valid.map(|parameters| (parameters.ungapped.clone(), gapped)),
        search_space: valid.map_or(0, |parameters| parameters.search_space),
        length_adjustment: valid.map_or(0, |parameters| i64::from(parameters.length_adjustment)),
        // An invalid context keeps `-searchsp`, else 0 (`effective_lengths`).
        context0_space: parameters[0].search_space,
    }
}

/// The defline line of the Bioseq: the line as read, else `>` and the title. BLASTX keeps
/// its own copy of `report::archive::defline_of`: its `FastaRecord` is its own type, and an
/// empty defline counts as missing.
fn defline_of(record: &FastaRecord) -> Cow<'_, [u8]> {
    match record
        .source
        .as_ref()
        .map(|source| source.defline.as_slice())
    {
        Some(defline) if !defline.is_empty() => Cow::Borrowed(defline),
        _ => {
            let mut line = vec![b'>'];
            line.extend_from_slice(record.title.as_bytes());
            Cow::Owned(line)
        }
    }
}

fn bioseq_input<'a>(record: &'a FastaRecord, defline: &'a [u8], protein: bool) -> BioseqInput<'a> {
    BioseqInput {
        id: &record.internal_id,
        defline,
        title: record.title.as_bytes(),
        letters: &record.sequence,
        protein,
    }
}

/// The frame of a context (`+1` .. `-3`).
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

/// The `masks` of the results: per query of a searched batch, one mask per frame with intervals,
/// in the order plus1, plus2, plus3, minus1, minus2, minus3; the intervals are the context's
/// nucleotide intervals other than an empty one and the whole query, in the order found (the
/// minus frames descending). The same intervals outfmt 15 and 16 print.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:930-959
/// (`Blast_GetSeqLocInfoVector`) and remote_blast.cpp:1064-1138.
fn results_masks(call: &ArchiveCall<'_>) -> Vec<MaskEntry> {
    let mut masks = Vec::new();
    for batch in call.batches {
        if batch.search_skipped {
            continue;
        }
        for qidx in 0..batch.queries.len() {
            let query = &call.queries[batch.query_ordinal + qidx];
            let length = query.sequence.len() as i32;
            for context in &batch.prepared.contexts[qidx * 6..qidx * 6 + 6] {
                let intervals: Vec<(i32, i32)> = context
                    .dna_masks
                    .iter()
                    .copied()
                    .filter(|&(from, to)| from <= to && !(from == 0 && to == length - 1))
                    .collect();
                if !intervals.is_empty() {
                    masks.push(MaskEntry {
                        id: query.internal_id.clone(),
                        frame: frame_of(context.frame),
                        intervals,
                    });
                }
            }
        }
    }
    masks
}

/// The `LCaseMask` parameters of `-lcase_masking`: per lower-case interval of a query (the
/// runs of the record, in nucleotide coordinates) one mask in frame `plus1` and one in `minus1`,
/// whatever `-strand`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:614-645
/// (the user mask of a translated query is mapped to its plus and minus strand).
fn lcase_masks(call: &ArchiveCall<'_>) -> Vec<MaskEntry> {
    let mut masks = Vec::new();
    for query in call.queries {
        for &(from, to) in &query.lowercase_masks {
            for frame in [Frame::Plus1, Frame::Minus1] {
                masks.push(MaskEntry {
                    id: query.internal_id.clone(),
                    frame,
                    intervals: vec![(from, to)],
                });
            }
        }
    }
    masks
}

/// The Seq-aligns of the batch: the elements of outfmt 8, one query at a time.
struct BlastxAligns<'a> {
    run: RunRecord<'a>,
    queries: &'a [FastaRecord],
    subjects: &'a [FastaRecord],
    names: &'a SeqalignSubjects<'a>,
    batches: &'a [BatchResults],
}

impl<'a> AlignSource for BlastxAligns<'a> {
    fn write<E: crate::report::serial::encoder::Encoder>(
        &mut self,
        out: &mut ArchiveOut<'_, E>,
    ) -> Result<()> {
        let mut segments = Segments::default();
        let mut subject_records: Vec<SubjectRecord<'a>> = Vec::new();
        let mut hsp_records: Vec<HspRecord<'a>> = Vec::new();
        for batch in self.batches {
            for (qidx, lists) in batch.queries.iter().enumerate() {
                let lists: &'a [HspList] = lists;
                if lists.iter().all(|list| list.hsps.is_empty()) {
                    continue;
                }
                let query = &self.queries[batch.query_ordinal + qidx];
                query_records(
                    lists,
                    self.subjects,
                    self.names,
                    &mut subject_records,
                    &mut hsp_records,
                )?;
                let query_record = QueryRecord {
                    id: SeqId::Local(&query.internal_id),
                    searched_length: i32::try_from(query.sequence.len())?,
                    offset: 0,
                };
                out.write_with(|enc, flush| {
                    write_seq_aligns(
                        enc,
                        &self.run,
                        &query_record,
                        &subject_records,
                        &hsp_records,
                        &mut segments,
                        |enc| flush(enc),
                    );
                });
            }
        }
        Ok(())
    }
}

/// Writes the archive of the batch of `call` to `writer`: the messages saved since the last
/// archive (the reader's warnings of the batch, the option warnings), the request, the
/// alignments of every query, the masks and the statistics of the first query.
pub(crate) fn write_blastx_archive(call: &ArchiveCall<'_>, writer: &mut dyn Write) -> Result<()> {
    let options = call.options;
    let run = RunRecord {
        program: Program::Blastx,
        db_title: match &call.names {
            SeqalignSubjects::Records => None,
            SeqalignSubjects::Database { title, .. } => Some(title),
        },
        gapped: options.gapped,
    };
    // The gapped block is the table's for the matrix and gap costs (LOSAT's BLASTX searches
    // BLOSUM62 with 11 / 1); an ungapped search has none.
    let gapped = options.gapped.then(|| {
        crate::stats::tables::lookup_protein_params(&crate::config::ProteinScoringSpec {
            matrix: crate::config::ScoringMatrix::Blosum62,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
        })
    });
    // The Bioseq inputs borrow their deflines: keep them alive for the whole call.
    let subject_deflines: Vec<Cow<'_, [u8]>> = call
        .subject_records
        .map(|records| records.iter().map(defline_of).collect())
        .unwrap_or_default();
    let subjects_input: Vec<BioseqInput<'_>> = call
        .subject_records
        .map(|records| {
            records
                .iter()
                .zip(&subject_deflines)
                .map(|(record, defline)| bioseq_input(record, defline, true))
                .collect()
        })
        .unwrap_or_default();
    let deflines: Vec<Cow<'_, [u8]>> = call.queries.iter().map(defline_of).collect();
    let queries: Vec<BioseqInput<'_>> = call
        .queries
        .iter()
        .zip(&deflines)
        .map(|(record, defline)| bioseq_input(record, defline, false))
        .collect();
    let first = call
        .batches
        .first()
        .map(|batch| first_context(batch, gapped))
        .unwrap_or(FirstContext {
            karlin: None,
            search_space: 0,
            length_adjustment: 0,
            context0_space: 0,
        });
    let searched = call.batches.iter().any(|batch| !batch.search_skipped);
    let ka_blocks = match &first.karlin {
        Some((ungapped, gapped)) => [
            ka_block(ungapped, false),
            gapped
                .as_ref()
                .map_or(KaBlock::missing(true), |gapped| ka_block(gapped, true)),
        ],
        None => [KaBlock::missing(false), KaBlock::missing(true)],
    };
    let lcase = lcase_masks(call);
    let masks = results_masks(call);
    let batch = BatchArchive {
        queries: &queries,
        subjects: &subjects_input,
        required: None,
        lcase_masks: &lcase,
        masks: &masks,
        ka_blocks,
        search_space: first.search_space,
        length_adjustment: first.length_adjustment,
        searched_space: searched.then_some(first.context0_space),
    };
    let mut aligns = BlastxAligns {
        run,
        queries: call.queries,
        subjects: call.subjects,
        names: &call.names,
        batches: call.batches,
    };
    write_archive(call.handle, &batch, &mut aligns, writer)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{try_parse_from, Cli, Commands};
    use crate::report::archive::options::{algorithm_options, ParamValue};
    use crate::report::serial::real::double_to_string_posix;

    // The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote for blastx (S/of8/b1/out, the OF8b1
    // probes; `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that each
    // archive adds: the name, the options given besides the inputs, and the list as
    // `Name=variant:value`.
    const CASES: &[(&str, &[&str], &str)] = &[
        ("bad", &[], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("batch_fast", &["-task", "blastx-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("batch_strand", &["-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbs0", &["-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbs0u", &["-comp_based_stats", "0u"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbs1", &["-comp_based_stats", "1"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:1 SmithWatermanMode=boolean:False"),
        ("cbs1u", &["-comp_based_stats", "1u"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:1 SmithWatermanMode=boolean:False"),
        ("cbs2", &["-comp_based_stats", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbs2u", &["-comp_based_stats", "2u"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbs3", &["-comp_based_stats", "3"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:False"),
        ("cbs3sw", &["-comp_based_stats", "3", "-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:True"),
        ("cbs3u", &["-comp_based_stats", "3u"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:False"),
        ("cbsD", &["-comp_based_stats", "D"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbsF", &["-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbst", &["-comp_based_stats", "t"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("gc11", &["-query_gencode", "11"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:11 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("gc2", &["-query_gencode", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:2 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("gc4_fast", &["-query_gencode", "4", "-task", "blastx-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:4 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("gc4", &["-query_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:4 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("lc_mask_all", &["-lcase_masking"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("lc_mask_fast", &["-lcase_masking", "-task", "blastx-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("lc_mask_gc4", &["-lcase_masking", "-query_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:4 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("lc_mask_minus", &["-lcase_masking", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("lc_mask_plus", &["-lcase_masking", "-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("lc_mask_seg", &["-lcase_masking", "-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:True WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("lc_mask_segno", &["-lcase_masking", "-seg", "no"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("mts_warn", &["-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:2 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_besthit", &["-subject_besthit"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 UseSubjectBestHit=boolean:True WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_cull0", &["-culling_limit", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 Culling=integer:0 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_cull5", &["-culling_limit", "5"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 Culling=integer:5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_cull_mts", &["-culling_limit", "3", "-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 Culling=integer:3 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:2 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_dbsize", &["-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ev", &["-evalue", "0.01"], "EvalueThreshold=cutoff:e-value:0.01 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_evbig", &["-evalue", "1000"], "EvalueThreshold=cutoff:e-value:1000 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_m80", &["-task", "blastx-fast", "-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 GapOpeningCost=integer:10 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WindowSize=integer:25 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_thr", &["-task", "blastx-fast", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:15 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_ws", &["-task", "blastx-fast", "-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 WordThreshold=integer:20 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gapext", &["-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapExtensionCost=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gapopen", &["-gapopen", "10"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gaps", &["-gapopen", "9", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_hsps", &["-max_hsps", "1"], "EvalueThreshold=cutoff:e-value:10 MaxHspsPerSubject=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_intron0", &["-max_intron_length", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_intron", &["-max_intron_length", "100"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:100 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_m45", &["-matrix", "BLOSUM45"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM45 WordThreshold=integer:15 WindowSize=integer:60 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_m62lc", &["-matrix", "blosum62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:blosum62 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_m80", &["-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WordThreshold=integer:13 WindowSize=integer:25 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_m80thr", &["-matrix", "BLOSUM80", "-threshold", "13"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WordThreshold=integer:13 WindowSize=integer:25 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_many", &["-evalue", "1e-5", "-word_size", "3", "-threshold", "12", "-window_size", "30", "-gapopen", "10", "-gapextend", "1", "-max_hsps", "2", "-qcov_hsp_perc", "30", "-query_gencode", "4", "-seg", "12 2.2 2.5", "-strand", "plus"], "EvalueThreshold=cutoff:e-value:1e-05 GapOpeningCost=integer:10 GapExtensionCost=integer:1 WordSize=integer:3 QueryCovHspPerc=real:30 MaxHspsPerSubject=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:30 UngappedMode=boolean:False QueryGeneticCode=integer:4 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_mts1", &["-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:1 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_mts4", &["-max_target_seqs", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:4 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_naln", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:3 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ndesc", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ndescaln", &["-num_descriptions", "3", "-num_alignments", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:4 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_pam30", &["-matrix", "PAM30"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM30 WordThreshold=integer:17 WindowSize=integer:15 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_pam70", &["-matrix", "PAM70"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM70 WordThreshold=integer:15 WindowSize=integer:20 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_qcov", &["-qcov_hsp_perc", "50"], "EvalueThreshold=cutoff:e-value:10 QueryCovHspPerc=real:50 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_soft", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:True SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_softf", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_sumstats", &["-sum_stats", "true"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:True LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_sumstatsf", &["-sum_stats", "false"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_thr_frac", &["-threshold", "11.7"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:11 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_thr", &["-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:15 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ungapped_cbs0", &["-ungapped", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:True QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("opt_ungapped_cbsF", &["-ungapped", "-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:True QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("opt_win0", &["-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:0 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_win", &["-window_size", "30"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:30 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws2", &["-word_size", "2"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws3", &["-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws4", &["-word_size", "4"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:4 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws5", &["-word_size", "5"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:19 WordSize=integer:5 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws6", &["-word_size", "6"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:21 WordSize=integer:6 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws7", &["-word_size", "7"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:20 WordSize=integer:7 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_xdrop", &["-xdrop_gap", "10", "-xdrop_gap_final", "30"], "EvalueThreshold=cutoff:e-value:10 GapXDropoff=real:10 GapXDropoffFinal=real:30 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_xdropu", &["-xdrop_ungap", "20"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("seg_no", &["-seg", "no"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("seg_par", &["-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("seg_plus", &["-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("seg_yes", &["-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:True WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("strand_both", &["-strand", "both"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("sw", &["-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:40 UngappedMode=boolean:False QueryGeneticCode=integer:1 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:True"),
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

    fn parse(extra: &[&str]) -> Option<BlastxArgs> {
        let mut argv = vec!["losat", "blastx", "-query", "q.fa", "-subject", "s.fa"];
        argv.extend_from_slice(extra);
        argv.extend_from_slice(&["-outfmt", "11"]);
        let cli: Cli = try_parse_from(argv).ok()?;
        match cli.command {
            Commands::Blastx(args) => Some(args),
            _ => None,
        }
    }

    /// The inputs of a command line that LOSAT's `blastx` refuses (a task, matrix, word size,
    /// composition mode or `-use_sw_tback` it does not search with): the options model must still
    /// write the list NCBI wrote for it. The options LOSAT takes come from its arguments, the
    /// refused ones are set by hand.
    fn refused_inputs(extra: &[&str]) -> OptionInputs {
        const REFUSED: [&str; 8] = [
            "-task",
            "-matrix",
            "-word_size",
            "-gapopen",
            "-gapextend",
            "-comp_based_stats",
            "-xdrop_gap",
            "-xdrop_gap_final",
        ];
        let mut kept = Vec::new();
        let mut refused: Vec<(&str, &str)> = Vec::new();
        let mut index = 0;
        while index < extra.len() {
            let flag = extra[index];
            if REFUSED.contains(&flag) {
                refused.push((flag, extra.get(index + 1).copied().unwrap_or_default()));
                index += 2;
            } else if matches!(flag, "-use_sw_tback" | "-xdrop_ungap") {
                refused.push((flag, extra.get(index + 1).copied().unwrap_or_default()));
                index += if flag == "-use_sw_tback" { 1 } else { 2 };
            } else {
                kept.push(flag);
                index += 1;
            }
        }
        let mut inputs = option_inputs(&parse(&kept).expect("the options LOSAT takes parse"));
        for (flag, value) in refused {
            match flag {
                "-task" if value == "blastx-fast" => inputs.task = "blastx-fast".to_string(),
                "-matrix" => inputs.matrix = Some(value.to_string()),
                "-word_size" => inputs.word_size = value.parse().ok(),
                "-gapopen" => inputs.gapopen = value.parse().ok(),
                "-gapextend" => inputs.gapextend = value.parse().ok(),
                "-xdrop_gap" => inputs.xdrop_gap = value.parse().ok(),
                "-xdrop_gap_final" => inputs.xdrop_gap_final = value.parse().ok(),
                "-comp_based_stats" => {
                    inputs.comp_based_stats = value.trim_end_matches('u').parse().unwrap_or(2);
                }
                "-use_sw_tback" => inputs.use_sw_tback = true,
                _ => {}
            }
        }
        inputs
    }

    #[test]
    fn the_option_lists_equal_those_of_ncbis_probes() {
        assert!(CASES.len() >= 70);
        let mut refused = 0;
        for &(name, extra, expected) in CASES {
            let inputs = match parse(extra) {
                Some(args) => option_inputs(&args),
                None => {
                    refused += 1;
                    refused_inputs(extra)
                }
            };
            assert_eq!(
                render(&algorithm_options(&inputs)),
                expected,
                "{name} {extra:?}"
            );
        }
        // The cases LOSAT's command line takes are the model's check of `option_inputs`.
        assert!(refused > 10 && CASES.len() - refused > 50, "{refused}");
    }

    #[test]
    fn the_service_is_plain() {
        assert_eq!(
            crate::report::archive::options::service(Program::Blastx, "blastx-fast"),
            "plain"
        );
    }
}
