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
use crate::core::blast_util::ncbi_int4_from_double;
use crate::report::archive::bioseq::BioseqInput;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{Filter, OptionInputs};
use crate::report::archive::requests::{LcaseFrames, RequestPolicy, RequestRecord, RunRequests};
use crate::report::archive::{
    ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive, KaBlock,
};
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::serial::seqalign::{write_seq_aligns, Segments};
use crate::stats::tables::KarlinParams;

/// The recorded options of a BLASTX run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2). RPSTBLASTN, whose search is BLASTX's, gives its own
/// (`BlastxArgs::request_inputs`), made from these and its arguments.
pub(crate) fn option_inputs(args: &BlastxArgs) -> OptionInputs {
    if let Some(inputs) = &args.request_inputs {
        return OptionInputs::clone(inputs);
    }
    let mut inputs = OptionInputs::new(Program::Blastx, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    inputs.evalue = Some(args.evalue);
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.matrix = args.matrix.clone();
    // `-xdrop_gap` and `-xdrop_gap_final` as given (`SetGapXDropoff`, `SetGapXDropoffFinal`);
    // `-xdrop_ungap` has no entry (NCBI records none).
    inputs.xdrop_gap = args.xdrop_gap;
    inputs.xdrop_gap_final = args.xdrop_gap_final;
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps;
    inputs.sum_stats = args.sum_stats;
    // `-threshold` is read as a double and recorded as an `int` (`static_cast<int>`), converted
    // as the x86-64 build does (`INT_MIN` past `Int4`, as the lookup table takes it:
    // `ncbi_int4_from_double`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:884-892
    // ```c++
    // CBlastOptions::SetWordThreshold(double w)
    // {
    //     ...
    //     if (m_Remote) {
    //         m_Remote->SetValue(eBlastOpt_WordThreshold, static_cast<int>(w));
    //     }
    // }
    // ```
    inputs.threshold = args.threshold.map(ncbi_int4_from_double);
    inputs.window_size = args.window_size;
    inputs.ungapped = args.ungapped;
    inputs.soft_masking = Some(args.soft_masking);
    // An archive is written only for a `-seg` that `resolve` read.
    inputs.seg = args.seg_spec().ok().map(Filter::from);
    inputs.culling_limit = args.culling_limit;
    inputs.subject_besthit = args.subject_besthit;
    inputs.max_intron_length = args.max_intron_length;
    inputs.query_gencode = i32::from(args.query_gencode);
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments.map(|size| size as i32))
        .map_or(500, |size| size);
    inputs.comp_based_stats = i32::from(super::args::composition_mode(&args.comp_based_stats));
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with BLASTX's recorded options).
pub(crate) fn handle_for(
    args: &BlastxArgs,
    archive_requested: bool,
    options: &crate::report::archive::OptionsHandle,
) -> Option<ArchiveHandle> {
    if !archive_requested {
        return None;
    }
    // NCBI keeps the name without the white space at its ends (`NStr::TruncateSpaces`,
    // blast_args.cpp:2472-2476), and the request's `database` is that name.
    let database = args.db.as_deref().map(crate::seqdb::db_argument_name);
    Some(ArchiveHandle::for_program(
        &option_inputs(args),
        database,
        options,
    ))
}

pub(crate) use crate::report::archive::error_lines;

/// No NCBI counterpart: the output streams of a library call in which a later batch may still be
/// refused by a BLAST XML, XML2 or JSON document (`local::report_batches`), held in the order in
/// which they are written: each write keeps its bytes with the index of its stream (`writer`), so
/// that `release` writes them in that order (the reports of every format, as a stream that they
/// all go to reads them from NCBI) and `discard` drops them. Outside `hold` the bytes pass
/// through. The side files of an outfmt 13 or 14 document of the run follow it (`side_files`):
/// held, written before the streams and dropped as they are.
pub(crate) struct ArchiveHold<S: Write>(std::sync::Arc<std::sync::Mutex<HoldState<S>>>);

struct HoldState<S> {
    holding: bool,
    discarded: bool,
    sinks: Vec<S>,
    held: Vec<(usize, Vec<u8>)>,
    /// The side files of the run's outfmt 13 and 14 documents.
    sides: Vec<super::blast_xml::HeldSideFiles>,
}

/// A stream of an `ArchiveHold`.
pub(crate) struct HoldWriter<S: Write> {
    hold: ArchiveHold<S>,
    sink: usize,
}

impl<S: Write> Clone for ArchiveHold<S> {
    fn clone(&self) -> Self {
        ArchiveHold(self.0.clone())
    }
}

impl<S: Write> ArchiveHold<S> {
    pub(crate) fn new() -> Self {
        ArchiveHold(std::sync::Arc::new(std::sync::Mutex::new(HoldState {
            holding: false,
            discarded: false,
            sinks: Vec::new(),
            held: Vec::new(),
            sides: Vec::new(),
        })))
    }

    fn state(&self) -> std::sync::MutexGuard<'_, HoldState<S>> {
        self.0
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// A stream that writes to `sink`.
    pub(crate) fn writer(&self, sink: S) -> HoldWriter<S> {
        let mut state = self.state();
        state.sinks.push(sink);
        HoldWriter {
            hold: self.clone(),
            sink: state.sinks.len() - 1,
        }
    }

    /// The side files `files` of an outfmt 13 or 14 document of the run, which follow this hold
    /// (`blast_xml::HeldSideFiles`): the command line holds its document's side files the same way
    /// (`blast_xml::HeldOutput`).
    pub(crate) fn side_files(
        &self,
        files: crate::report::blast_xml2::RunSideFiles,
    ) -> crate::report::blast_xml2::RunSideFiles {
        let side = super::blast_xml::HeldSideFiles::new(files);
        let mut state = self.state();
        if state.discarded {
            side.discard();
        } else if state.holding {
            side.hold();
        }
        state.sides.push(side.clone());
        Box::new(side)
    }

    /// Keeps everything written from now on.
    pub(crate) fn hold(&self) {
        let mut state = self.state();
        state.holding = true;
        for side in &state.sides {
            side.hold();
        }
    }

    /// Writes what is held and passes the bytes on again: the side files first (in order), then
    /// the streams in the order of their writes. A side file that cannot be opened stops the run
    /// there, as NCBI's does (`Cannot open output file`, exit status 1): the side files before it
    /// stay, and the streams and the side files of the later documents get nothing of what was
    /// held.
    pub(crate) fn release(&self) -> Result<()> {
        let mut state = self.state();
        state.holding = false;
        let held = std::mem::take(&mut state.held);
        if state.discarded {
            return Ok(());
        }
        for (index, side) in state.sides.iter().enumerate() {
            if let Err(error) = side.release() {
                // Nothing more of what was held is written: neither the streams, nor the side
                // files of the documents after this one.
                for later in &state.sides[index + 1..] {
                    later.discard();
                }
                state.discarded = true;
                return Err(error);
            }
        }
        for (sink, bytes) in held {
            state.sinks[sink]
                .write_all(&bytes)
                .map_err(super::native::output_error)?;
        }
        for sink in &mut state.sinks {
            sink.flush().map_err(super::native::output_error)?;
        }
        Ok(())
    }

    /// Drops what is held and everything written after.
    pub(crate) fn discard(&self) {
        let mut state = self.state();
        state.discarded = true;
        state.held.clear();
        for side in &state.sides {
            side.discard();
        }
    }

    /// The end of a run that stopped with `error` and held its output for its BLAST XML, XML2 or
    /// JSON documents: dropped after one of the refusals that such a hold is for
    /// (`blast_xml::drops_held_document`), written after any other error (the bytes of the batches
    /// before it, as the run without a hold writes them), as the command line's document
    /// (`blast_xml::HeldOutput`).
    pub(crate) fn end_document_with(&self, error: &anyhow::Error) {
        if super::blast_xml::drops_held_document(error) {
            self.discard();
        } else {
            let _ = self.release();
        }
    }
}

impl<S: Write> Write for HoldWriter<S> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let mut state = self.hold.state();
        if state.discarded {
            return Ok(bytes.len());
        }
        if state.holding {
            match state.held.last_mut() {
                Some((sink, held)) if *sink == self.sink => held.extend_from_slice(bytes),
                _ => state.held.push((self.sink, bytes.to_vec())),
            }
            return Ok(bytes.len());
        }
        state.sinks[self.sink].write(bytes)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let mut state = self.hold.state();
        if state.discarded || state.holding {
            return Ok(());
        }
        state.sinks[self.sink].flush()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
// ```c++
// CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
//     : m_QueryVector(& queries)
// {
//     if (queries.Empty()) {
//         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:650-652
// ```c++
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
/// Whether NCBI searches the reader batch `batch` and writes its archive: a batch without a
/// query (every record skipped by `-query_loc`) is `Empty CBlastQueryVector`, and one without a
/// valid context and with a query without letters fails its set-up
/// (`ReaderBatchShape::set_up_fails`), which stops the run there as the search does
/// (`native::run_inner`, `query_setup::prepare_unit_queries`).
pub(crate) fn written_batch(batch: &[FastaRecord], options: &ResolvedOptions) -> bool {
    !batch.is_empty()
        && !super::reader_batch::ReaderBatchShape::layout(batch, options).set_up_fails()
}

/// The queries of one query batch and what its archive reads (`with_batch`).
pub(crate) struct CallBatch<'a> {
    pub options: &'a ResolvedOptions,
    /// The batch's queries; `BatchResults::query_ordinal` indexes them.
    pub queries: &'a [FastaRecord],
    /// The requests of the run's archives (`requests`, with every `-subject` record; made once).
    pub requests: &'a RunRequests<'a>,
    /// The results of the batch: one entry, or several when the search split a batch above
    /// 10002 letters (`BATCH_SIZE`).
    pub batches: &'a [BatchResults],
}

/// What the archive of one query batch needs.
pub(crate) struct ArchiveCall<'a> {
    pub handle: &'a ArchiveHandle,
    pub batch: CallBatch<'a>,
    /// The subjects by the oids of the hit lists.
    pub subjects: &'a [FastaRecord],
    pub names: SeqalignSubjects<'a>,
}

/// The values of a batch's first query that `search-stats` and `ka-blocks` take (the
/// `EffectiveSearchSpace` option: `batch_searched_space`).
struct FirstContext {
    /// The ungapped block of the first valid context and the gapped one (absent for an
    /// ungapped search); `None` when the first query has no valid context (both blocks are
    /// then missing).
    karlin: Option<(KarlinParams, Option<KarlinParams>)>,
    /// The effective search space and length adjustment of that context (0 without one).
    search_space: i64,
    length_adjustment: i64,
}

fn first_context(batch: &BatchResults, gapped: Option<KarlinParams>) -> FirstContext {
    let parameters = &batch.parameters[..6];
    let valid = parameters.iter().find(|parameters| parameters.valid);
    FirstContext {
        karlin: valid.map(|parameters| (parameters.ungapped, gapped)),
        search_space: valid.map_or(0, |parameters| parameters.search_space),
        length_adjustment: valid.map_or(0, |parameters| i64::from(parameters.length_adjustment)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1718-1725
// ```c++
// inline Int8
// CBlastOptionsLocal::GetEffectiveSearchSpace() const
// {
//     if (m_EffLenOpts->num_searchspaces == 0)
//         return 0;
//
//     return m_EffLenOpts->searchsp_eff[0];
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
// ```c++
//     int status = m_PrelimSearch->CheckInternalData();
//     if (status != 0)
//     {
//          // Search was not run, but we send back an empty CSearchResultSet.
// ```
/// The `EffectiveSearchSpace` option of a reader batch whose results are `batches` (its units;
/// `BatchArchive::searched_space`, `StrategyExport::note_searched`): the search space of the first
/// query's context 0 (frame +1; an invalid context keeps `-searchsp`, else 0,
/// `effective_lengths`), when NCBI searched the batch (a unit was searched).
pub(crate) fn batch_searched_space(batches: &[BatchResults]) -> Option<i64> {
    batches
        .iter()
        .any(|batch| !batch.search_skipped)
        .then(|| batches[0].parameters[0].search_space)
}

/// The requests of a run's batches (`RunRequests`): nucleotide queries, protein subjects (the
/// `-subject` records `subjects`, `None` for `-db`), and the lower-case runs of each query in both
/// frames of a translated query (the reader keeps them only with `-lcase_masking`,
/// `FastaRecord::archive_lowercase_runs`).
pub(crate) fn requests(subjects: Option<&[FastaRecord]>) -> RunRequests<'_> {
    RunRequests::new(
        RequestPolicy {
            query_protein: false,
            subject_protein: true,
            lcase: Some(LcaseFrames::PlusAndMinus),
        },
        subjects,
    )
}

impl RequestRecord for FastaRecord {
    fn bioseq(&self, protein: bool) -> BioseqInput<'_> {
        // The defline line: the line as read, else `>` and the title (`report::archive::defline_of`;
        // BLASTX's `FastaRecord` is its own type, and an empty defline counts as missing).
        let defline = match self.source.as_ref().map(|source| source.defline.as_slice()) {
            Some(defline) if !defline.is_empty() => Cow::Borrowed(defline),
            _ => {
                let mut line = vec![b'>'];
                line.extend_from_slice(self.title.as_bytes());
                Cow::Owned(line)
            }
        };
        BioseqInput {
            id: &self.internal_id,
            ids: self.believed().map(|ids| ids.ids.as_slice()),
            seq_loc: self.archive_seq_loc(),
            defline,
            title: self.title.as_bytes(),
            // The whole record, also when a range cut it (`FastaRecord::archive_letters`).
            letters: self.archive_letters(),
            protein,
        }
    }

    fn searched_len(&self) -> usize {
        self.sequence.len()
    }

    /// The runs of the whole record, also outside a `-query_loc` interval (the reader's mask of
    /// the record, `x_FastaToSeqLoc`), whatever `-strand`.
    fn lowercase_runs(&self) -> Vec<(i32, i32)> {
        self.archive_lowercase_runs()
    }

    fn mask(&self, frame: Frame, intervals: Vec<(i32, i32)>) -> MaskEntry {
        MaskEntry {
            id: self.internal_id.clone(),
            best: self.believed().map(|ids| ids.best.clone()),
            frame,
            intervals,
        }
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
/// nucleotide intervals moved into the record (`seq_range::record_mask`: other than an empty one
/// and the whole interval), in the order found (the minus frames descending). The same intervals
/// outfmt 15 and 16 print.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:930-959
/// (`Blast_GetSeqLocInfoVector`) and remote_blast.cpp:1064-1138.
fn results_masks(call: &CallBatch<'_>) -> Vec<MaskEntry> {
    let mut masks = Vec::new();
    for batch in call.batches {
        if batch.search_skipped {
            continue;
        }
        for qidx in 0..batch.queries.len() {
            let query = &call.queries[batch.query_ordinal + qidx];
            let (length, offset) = (query.sequence.len(), query.placement().offset);
            for context in &batch.prepared.contexts[qidx * 6..qidx * 6 + 6] {
                let intervals: Vec<(i32, i32)> = context
                    .dna_masks
                    .iter()
                    .filter_map(|&(from, to)| {
                        crate::blastinput::seq_range::record_mask(from, to, length, offset)
                    })
                    .collect();
                if !intervals.is_empty() {
                    masks.push(MaskEntry {
                        id: query.internal_id.clone(),
                        best: query.believed().map(|ids| ids.best.clone()),
                        frame: frame_of(context.frame),
                        intervals,
                    });
                }
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
                    id: query.seq_id(),
                    searched_length: i32::try_from(query.sequence.len())?,
                    offset: i32::try_from(query.placement().offset)?,
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

/// Calls `each` with the archive of the batch of `call` (`BatchArchive`: its request,
/// `RunRequests::with_batch`, and its results).
fn with_batch(
    call: &CallBatch<'_>,
    each: impl FnOnce(&BatchArchive<'_>) -> Result<()>,
) -> Result<()> {
    let options = call.options;
    // The gapped block is the table's for the matrix and gap costs (LOSAT's BLASTX searches
    // BLOSUM62 with 11 / 1); an ungapped search has none.
    let gapped = options.gapped.then(|| {
        crate::stats::tables::lookup_protein_params(&crate::config::ProteinScoringSpec {
            matrix: options.scoring,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
        })
    });
    let first = call
        .batches
        .first()
        .map(|batch| first_context(batch, gapped))
        .unwrap_or(FirstContext {
            karlin: None,
            search_space: 0,
            length_adjustment: 0,
        });
    let ka_blocks = match &first.karlin {
        Some((ungapped, gapped)) => [
            ka_block(ungapped, false),
            gapped
                .as_ref()
                .map_or(KaBlock::missing(true), |gapped| ka_block(gapped, true)),
        ],
        None => [KaBlock::missing(false), KaBlock::missing(true)],
    };
    let masks = results_masks(call);
    call.requests.with_batch(
        call.queries,
        call.queries.first().map(FastaRecord::placement),
        |request| {
            each(&BatchArchive {
                request: *request,
                masks: &masks,
                ka_blocks,
                search_space: first.search_space,
                length_adjustment: first.length_adjustment,
                searched_space: batch_searched_space(call.batches),
            })
        },
    )
}

/// Writes the archive of the batch of `call` to `writer`: the messages saved since the last
/// archive (the reader's warnings of the batch, the option warnings), the request, the
/// alignments of every query, the masks and the statistics of the first query. NCBI's posts of
/// the bytes that the archive's stream replaces go to `diagnostics` (`write_archive`).
pub(crate) fn write_blastx_archive(
    call: &ArchiveCall<'_>,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Blastx,
        db_title: match &call.names {
            SeqalignSubjects::Records => None,
            SeqalignSubjects::Database { title, .. } => Some(title),
        },
        gapped: call.batch.options.gapped,
    };
    with_batch(&call.batch, |batch| {
        let mut aligns = BlastxAligns {
            run,
            queries: call.batch.queries,
            subjects: call.subjects,
            names: &call.names,
            batches: call.batch.batches,
        };
        write_archive(call.handle, batch, &mut aligns, writer, diagnostics)
    })
}

#[cfg(test)]
mod hold_tests {
    //! `ArchiveHold`: what a run that held its output writes when it ends (`end_document_with`:
    //! the XML, XML2 and JSON documents), and the side files of outfmt 13 and 14 that follow it
    //! (`side_files`).
    use super::*;
    use crate::cli::NativeError;
    use crate::report::blast_xml2::SideFiles;
    use crate::seqdb::SeqDbError;
    use std::sync::{Arc, Mutex};

    /// A stream that the test reads after the hold wrote to it.
    #[derive(Clone, Default)]
    struct Shared(Arc<Mutex<Vec<u8>>>);

    impl Write for Shared {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(bytes.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    impl Shared {
        fn text(&self) -> String {
            String::from_utf8(self.0.lock().unwrap().clone()).unwrap()
        }
    }

    /// Side files that record each file as `<n>:<bytes>;`.
    struct Recorded(Shared);

    impl SideFiles for Recorded {
        fn create(&mut self, _n: u32) -> std::io::Result<Box<dyn Write + '_>> {
            Err(std::io::Error::other("a run writes through write_file"))
        }

        fn listed(&self, n: u32) -> Vec<u8> {
            format!("f{n}").into_bytes()
        }

        fn write_file(&mut self, n: u32, report: &[u8]) -> Result<()> {
            write!(self.0, "{n}:{};", String::from_utf8_lossy(report))?;
            Ok(())
        }
    }

    /// A run that holds two streams and its side files, wrote `a1 b1 a2` and two side files, and
    /// has not ended.
    fn held_run() -> (ArchiveHold<Shared>, [Shared; 3]) {
        let (a, b, files) = (Shared::default(), Shared::default(), Shared::default());
        let hold = ArchiveHold::new();
        let (mut first, mut second) = (hold.writer(a.clone()), hold.writer(b.clone()));
        hold.hold();
        let mut side = hold.side_files(Box::new(Recorded(files.clone())));
        first.write_all(b"a1 ").unwrap();
        side.write_file(1, b"one").unwrap();
        second.write_all(b"b1 ").unwrap();
        first.write_all(b"a2 ").unwrap();
        side.write_file(2, b"two").unwrap();
        assert_eq!(side.listed(2), b"f2");
        for stream in [&a, &b, &files] {
            assert_eq!(stream.text(), "", "nothing passes a hold");
        }
        (hold, [a, b, files])
    }

    fn written(streams: &[Shared; 3]) -> [String; 3] {
        [streams[0].text(), streams[1].text(), streams[2].text()]
    }

    const ALL: [&str; 3] = ["a1 a2 ", "b1 ", "1:one;2:two;"];
    const NOTHING: [&str; 3] = ["", "", ""];

    fn ncbis_error() -> anyhow::Error {
        NativeError {
            exit: 3,
            message: "BLAST engine error: Warning: Sequence contains no data\n".into(),
        }
        .into()
    }

    /// A refusal of the taxonomy database as the XML2 and JSON documents meet it (`XmlRefusal`,
    /// made of the reader's `TaxonomyRefusal`, which outfmt 5's prefetch returns as it is).
    fn document_refusal() -> anyhow::Error {
        super::super::blast_xml::XmlRefusal::of_taxonomy(
            crate::algorithm::common::database::taxdb_error_for_tests(SeqDbError::Unsupported {
                path: "taxdb.btd".into(),
                message: "the taxonomy record of taxid 9606 reaches past the mapped end".into(),
            }),
        )
    }

    #[test]
    fn a_release_writes_the_side_files_and_every_stream_in_order() {
        let (hold, streams) = held_run();
        hold.release().unwrap();
        assert_eq!(written(&streams), ALL);
    }

    /// NCBI's error of a side file that it cannot open (`CArgException`, exit status 1): no
    /// error archive follows it (`error_lines`), and NCBI wrote the bytes before it all the same.
    fn side_file_error() -> anyhow::Error {
        crate::report::blast_xml2::framing::side_file_open_error()
    }

    #[test]
    fn the_documents_hold_drops_after_their_refusal_and_writes_after_any_other_error() {
        let (hold, streams) = held_run();
        hold.end_document_with(&document_refusal());
        assert_eq!(written(&streams), NOTHING);
        // The taxonomy database's refusal of the prefetch (outfmt 5) is one of them.
        let (hold, streams) = held_run();
        let prefetch =
            crate::algorithm::common::database::taxdb_error_for_tests(SeqDbError::Unsupported {
                path: "taxdb.btd".into(),
                message: "the taxonomy record of taxid 9606 reaches past the mapped end".into(),
            });
        hold.end_document_with(&prefetch);
        assert_eq!(written(&streams), NOTHING);
        for other in [
            ncbis_error(),
            side_file_error(),
            anyhow::anyhow!("taxdb.btd: Input/output error; this database is not supported"),
        ] {
            let (hold, streams) = held_run();
            hold.end_document_with(&other);
            assert_eq!(written(&streams), ALL, "{other}");
        }
    }

    /// Side files whose every file cannot be opened.
    struct Blocked;

    impl SideFiles for Blocked {
        fn create(&mut self, _n: u32) -> std::io::Result<Box<dyn Write + '_>> {
            Err(std::io::Error::other("blocked"))
        }

        fn listed(&self, _n: u32) -> Vec<u8> {
            Vec::new()
        }
    }

    /// The side files of one document cannot be opened at the release (they could when their
    /// reports were made): NCBI's error ends the run, and neither the streams nor the side files
    /// of the document after it get what was held, then or when the run's writers are dropped.
    #[test]
    fn a_side_file_error_at_the_release_writes_nothing_after_it() {
        let (a, files) = (Shared::default(), Shared::default());
        let hold = ArchiveHold::new();
        let mut stream = hold.writer(a.clone());
        hold.hold();
        let mut first = hold.side_files(Box::new(Recorded(files.clone())));
        let mut blocked = hold.side_files(Box::new(Blocked));
        let mut last = hold.side_files(Box::new(Recorded(files.clone())));
        stream.write_all(b"a1 ").unwrap();
        first.write_file(1, b"one").unwrap();
        blocked.write_file(1, b"x").unwrap();
        last.write_file(1, b"late").unwrap();
        let error = hold.release().unwrap_err();
        let native = error.downcast_ref::<NativeError>().expect("NCBI's error");
        assert_eq!(
            native.message,
            "Command line argument error: Cannot open output file\n"
        );
        drop((first, blocked, last, stream, hold));
        assert_eq!([a.text(), files.text()], ["", "1:one;"].map(String::from));
    }

    #[test]
    fn after_a_discard_nothing_more_is_written() {
        let (hold, streams) = held_run();
        let mut late = hold.writer(streams[0].clone());
        let mut side = hold.side_files(Box::new(Recorded(streams[2].clone())));
        hold.discard();
        late.write_all(b"late").unwrap();
        side.write_file(3, b"three").unwrap();
        hold.release().unwrap();
        assert_eq!(written(&streams), NOTHING);
    }
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
                    ParamValue::IntegerList(values) => format!("integer-list:{values:?}"),
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

    /// The inputs of a command line that LOSAT's `blastx` parser refuses (`-use_sw_tback`, which it
    /// does not search with): the options model must still write the list NCBI
    /// wrote for it. The options LOSAT takes come from its arguments, the refused ones are set by
    /// hand.
    fn refused_inputs(extra: &[&str]) -> OptionInputs {
        const REFUSED: [&str; 6] = [
            "-task",
            "-matrix",
            "-word_size",
            "-gapopen",
            "-gapextend",
            "-comp_based_stats",
        ];
        let mut kept = Vec::new();
        let mut refused: Vec<(&str, &str)> = Vec::new();
        let mut index = 0;
        while index < extra.len() {
            let flag = extra[index];
            if REFUSED.contains(&flag) {
                refused.push((flag, extra.get(index + 1).copied().unwrap_or_default()));
                index += 2;
            } else if flag == "-use_sw_tback" {
                refused.push((flag, ""));
                index += 1;
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
                "-comp_based_stats" => {
                    inputs.comp_based_stats =
                        i32::from(crate::algorithm::blastx::args::composition_mode(value));
                }
                "-use_sw_tback" => inputs.use_sw_tback = true,
                _ => {}
            }
        }
        inputs
    }

    // A `-threshold` past `Int4` (2147483648, 1e10, +inf, 1e999) is recorded as `INT_MIN`, as
    // NCBI BLAST+ 2.17.0 converts it (the lookup table takes the same value). NCBI's build writes
    // it as `-8` in the archive, from a signed overflow that LOSAT does not reproduce
    // (`report/serial/asn_text.rs` `itoa_i32`).
    #[test]
    fn a_threshold_past_int4_is_recorded_as_int_min() {
        let inputs = |value: &str| {
            option_inputs(&parse(&["-threshold", value]).unwrap_or_else(|| panic!("{value}")))
        };
        let listed = |threshold: &str| {
            let case = CASES
                .iter()
                .find(|case| case.0 == "opt_thr_frac")
                .expect("opt_thr_frac");
            case.2.replace(
                "WordThreshold=integer:11",
                &format!("WordThreshold=integer:{threshold}"),
            )
        };
        for value in ["2147483648", "1e10", "+inf", "1e999"] {
            assert_eq!(
                render(&algorithm_options(&inputs(value))),
                listed("-2147483648"),
                "{value}"
            );
        }
        assert_eq!(
            render(&algorithm_options(&inputs("2147483647.5"))),
            listed("2147483647")
        );
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
        // The cases LOSAT's command line takes are the model's check of `option_inputs` (every
        // matrix and composition mode since series SC, SC6; the word sizes LOSAT does not search
        // with since line LIB, P-S3, whose parser reads them and whose `resolve` refuses them):
        // 2 of the 80 are refused (`-use_sw_tback`; the X-drop cases parse since line RG X1).
        assert!(refused == 2 && CASES.len() - refused > 60, "{refused}");
    }

    #[test]
    fn the_service_is_plain() {
        assert_eq!(
            crate::report::archive::options::service(Program::Blastx, "blastx-fast"),
            "plain"
        );
    }
}
