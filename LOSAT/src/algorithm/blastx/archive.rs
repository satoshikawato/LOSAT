//! BLASTX's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:283-286
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
    // The archive's `HitlistSize` is that of outfmt 11, whoever made the other options.
    let hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments.map(|size| size as i32))
        .map_or(args.formatting.recorded_hitlist_start(), |size| size);
    if let Some(inputs) = &args.request_inputs {
        let mut inputs = OptionInputs::clone(inputs);
        inputs.hitlist_size = hitlist_size;
        return inputs;
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
    inputs.imported = args.import.as_ref().map(|run| run.recorded.clone());
    inputs.hitlist_size = hitlist_size;
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:650-652
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
    // RPSTBLASTN's first query of a batch has its profile's K (`BatchResults::gapped_k`, the
    // footer's and the XML's), the table's otherwise.
    let gapped = gapped.map(|mut block| {
        if let Some(&k) = batch
            .gapped_k
            .as_ref()
            .and_then(|gapped_k| gapped_k.first())
        {
            block.k = k;
        }
        block
    });
    FirstContext {
        karlin: valid.map(|parameters| (parameters.ungapped, gapped)),
        search_space: valid.map_or(0, |parameters| parameters.search_space),
        length_adjustment: valid.map_or(0, |parameters| i64::from(parameters.length_adjustment)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1718-1725
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
/// The `EffectiveSearchSpace` option of a reader batch whose results are `batches` (its units;
/// `BatchArchive::searched_space`, `StrategyExport::note_searched`): the search space of the first
/// query's context 0 (frame +1; an invalid context keeps `-searchsp`, else 0,
/// `effective_lengths`; before the batch's first context, `-strand minus`, the option is
/// `-searchsp` too: NCBI's `SplitQuery_SetEffectiveSearchSpace` leaves a set option as it is)
/// and that context's query length (0 for an invalid one; `archive::RequestLength::space`), when
/// NCBI searched the batch (a unit was searched); for RPSTBLASTN the space its last search unit
/// leaves in the options handle (`BatchResults::request_space`, which holds `-searchsp` as given
/// and 0 for an invalid context 0).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:157-161
pub(crate) fn batch_searched_space(
    batches: &[BatchResults],
    options: &ResolvedOptions,
) -> Option<crate::report::archive::SearchedContext0> {
    batches.iter().any(|batch| !batch.search_skipped).then(|| {
        let parameters = &batches[0].parameters[0];
        crate::report::archive::SearchedContext0 {
            space: batches[0]
                .request_space
                .unwrap_or(if batches[0].prepared.first_context > 0 {
                    options.effective_search_space
                } else {
                    parameters.search_space
                }),
            query_length: if parameters.valid {
                batches[0].prepared.contexts[0].length as i64
            } else {
                0
            },
        }
    })
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
                line.extend_from_slice(&self.title);
                Cow::Owned(line)
            }
        };
        BioseqInput {
            id: &self.internal_id,
            ids: self.believed().map(|ids| ids.ids.as_slice()),
            seq_loc: self.archive_seq_loc(),
            defline,
            title: &self.title,
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
                searched_space: batch_searched_space(call.batches, options),
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
    // `-subject`'s Seq-loc mapper cuts an HSP of length 0 at the first letter of the searched part
    // of a subject with HSPs of letters in two, and makes alignments without scores of a subject
    // whose every HSP has length 0 where it rewrites one, which LOSAT does not write
    // (`empty_rows::refuse_mapped_without_scores`).
    if run.remaps_subjects() && !run.gapped {
        super::empty_rows::refuse_mapped_without_scores(call.batch.batches, run.gapped, None)?;
    }
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

/// The arguments a blastx run keeps for another program's strategy (`import::RunAs`).
fn carried(args: BlastxArgs) -> crate::report::archive::read::import::Carried {
    crate::report::archive::read::import::Carried {
        app: "blastx",
        path: args.import_search_strategy.expect("an import"),
        source: args.import.as_ref().and_then(|run| run.source.clone()),
        query: args.query,
        query_loc: args.query_loc,
        strand: args.strand,
        lcase_masking: args.lcase_masking,
        parse_deflines: args.parse_deflines,
        id_scope: None,
        out: args.out,
        outfmt: args.outfmt,
        show_gis: args.show_gis,
        num_descriptions: args.num_descriptions,
        num_alignments: args.num_alignments,
        max_target_seqs: args.max_target_seqs.map(i64::from),
        formatting: args.formatting,
        num_threads: usize::try_from(args.num_threads).unwrap_or(1),
        mt_mode: args.mt_mode,
        db: args.db,
        subject: args.subject,
        subject_loc: args.subject_loc,
        db_filter: args.db_filter,
        dbsize: args.dbsize,
    }
}

/// blastx's arguments for a BLASTX strategy that another program runs (`import::RunAs`): its
/// defaults with the arguments the running program keeps.
pub(crate) fn from_carried(
    carried: crate::report::archive::read::import::Carried,
) -> Result<BlastxArgs> {
    let cli: crate::cli::Cli = crate::cli::try_parse_from(["losat", "blastx"])?;
    let crate::cli::Commands::Blastx(mut out) = cli.command else {
        unreachable!("the arguments of blastx")
    };
    out.import = carried.import_app();
    out.max_target_seqs = carried
        .max_target_seqs
        .map(|n| i32::try_from(n).unwrap_or(i32::MAX));
    out.num_threads = i32::try_from(carried.num_threads).unwrap_or(i32::MAX);
    out.import_search_strategy = Some(carried.path);
    out.query = carried.query;
    out.query_loc = carried.query_loc;
    out.strand = carried.strand;
    out.lcase_masking = carried.lcase_masking;
    out.parse_deflines = carried.parse_deflines;
    out.out = carried.out;
    out.outfmt = carried.outfmt;
    out.show_gis = carried.show_gis;
    out.num_descriptions = carried.num_descriptions;
    out.num_alignments = carried.num_alignments;
    out.formatting = carried.formatting;
    out.mt_mode = carried.mt_mode;
    out.db = carried.db;
    out.subject = carried.subject;
    out.subject_loc = carried.subject_loc;
    out.db_filter = carried.db_filter;
    out.dbsize = carried.dbsize;
    Ok(out)
}

/// `-import_search_strategy` of BLASTX (`RecoverSearchStrategy`, `SetOptionsForSavedStrategy`,
/// blast_args.cpp:3800-3836, quoted at `blastp::archive::apply_strategy`): the arguments of a
/// command line that gives nothing (the strategy task's defaults) with the strategy's options set
/// on them in the file's order, and from `args` only what NCBI's re-read argument classes take:
/// the query options, `-out`, the formatting options, the threads and the database arguments
/// when `-db` or `-subject` is given. The run then goes on as any other (one search path).
pub(crate) fn apply_strategy(
    args: BlastxArgs,
    diagnostics: &mut dyn std::io::Write,
) -> Result<BlastxArgs> {
    use crate::report::archive::read::import::{self, Command, Start};
    use crate::report::archive::read::options::{EProgram, Setter};
    let program = Program::Blastx.name();
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
        args.import.as_ref().map_or("blastx", |run| run.app),
        diagnostics,
    )?;
    if start.imported.program() != EProgram::Blastx {
        return Err(start.run_as(program, || carried(args)));
    }
    let cli: crate::cli::Cli = crate::cli::try_parse_from(["losat", "blastx"])?;
    let crate::cli::Commands::Blastx(mut out) = cli.command else {
        unreachable!("the arguments of blastx")
    };
    let mut run = start.run(args.parse_deflines, program)?;
    run.query_range = start.query_range(args.query_loc.as_ref(), program)?;
    // `CQueryOptionsArgs`, `CStdCmdLineArgs`, `CFormattingArgs`, `CMTArgs`.
    out.query = args.query;
    out.query_loc = args.query_loc;
    out.strand = args.strand;
    out.lcase_masking = args.lcase_masking;
    out.parse_deflines = args.parse_deflines;
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
            Setter::EvalueThreshold(v) => out.evalue = *v,
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
            Setter::MaxHspsPerSubject(n) => out.max_hsps = Some(*n).filter(|n| *n > 0),
            Setter::QueryCovHspPerc(v) => out.qcov_hsp_perc = Some(*v),
            Setter::MaskAtHash(on) => out.soft_masking = *on,
            Setter::GappedMode(gapped) => out.ungapped = !*gapped,
            Setter::SumStatisticsMode(on) => out.sum_stats = Some(*on),
            Setter::QueryGeneticCode(n) => {
                out.query_gencode = u8::try_from(*n).map_err(|_| {
                    Start::unported("a search strategy with this query genetic code", program)
                })?
            }
            Setter::LongestIntronLength(n) => out.max_intron_length = *n,
            Setter::SubjectBestHit => out.subject_besthit = true,
            Setter::DustFiltering(false)
            | Setter::RepeatFiltering(false)
            | Setter::SmithWatermanMode(false)
            | Setter::OutOfFrameMode(false) => {}
            other => return Err(Start::unported_setter(other, program)),
        }
    }
    if let Some(value) = seg.value() {
        out.seg = value;
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
