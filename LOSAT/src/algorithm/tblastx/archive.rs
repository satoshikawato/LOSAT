//! TBLASTX's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:198
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

use std::io::Write;
use std::ops::Range;

use anyhow::Result;

use super::blast_engine::{TblastxHsp, TblastxQueryStats};
use super::report::SubjectLetters;
use super::seqalign::{query_records_of, SeqalignSubjects};
use super::subjects::SubjectSetNumbers;
use super::TblastxArgs;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::core::blast_util::ncbi_int4_from_double;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{Filter, OptionInputs};
use crate::report::archive::read::import;
use crate::report::archive::read::options::{EProgram, Setter};
use crate::report::archive::requests::{LcaseFrames, RequestPolicy, RunRequests};
use crate::report::archive::{
    ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive, KaBlock,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::serial::seqalign::{write_seq_aligns, Segments};

/// What a TBLASTX run's arguments give its search strategy's requests
/// (`strategy::StrategyExport::for_args`): the archive's recorded options, the counts
/// of the `HitlistSize` and the database of `-db`.
impl crate::report::archive::strategy::StrategyArgs for TblastxArgs {
    fn strategy_inputs(&self) -> crate::report::archive::strategy::StrategyInputs<'_> {
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

/// The recorded options of a TBLASTX run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2). TBLASTX has one task, and the options it does not
/// accept (`-entrez_query`, ...) are never recorded. `-strand` is not an option of the request:
/// NCBI keeps it in the queries' Seq-locs (the query reader's strand, `with_strand`).
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
    // an `int` (`static_cast<int>`), converted as the x86-64 build does (`INT_MIN` past `Int4`,
    // as the lookup table takes it: `ncbi_int4_from_double`). Without `-threshold` and
    // `-window_size` the list holds the matrix's suggestions (`report::archive::options`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:884-892
    inputs.matrix = args.matrix.clone();
    inputs.threshold = args.threshold.map(ncbi_int4_from_double);
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
    // `MaskAtHash` is recorded with the value given (the class always sets it: false when absent).
    inputs.soft_masking = args.soft_masking;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:321-323
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:2043-2052
    // `SumStatistics` is recorded only when `-sum_stats` is given, with its value.
    inputs.sum_stats = args.sum_stats;
    inputs.window_size = args.window_size;
    inputs.seg = args.seg_spec().ok().map(Filter::from);
    inputs.culling_limit = args.culling_limit;
    inputs.best_hit_overhang = args.best_hit_overhang;
    inputs.best_hit_score_edge = args.best_hit_score_edge;
    inputs.subject_besthit = args.subject_besthit;
    inputs.query_gencode = i32::from(args.query_gencode);
    inputs.db_gencode = i32::from(args.db_gencode);
    inputs.imported = args.import.as_ref().map(|run| run.recorded.clone());
    inputs.application = args.import.as_ref().map(|run| run.app.ncbi_name());
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:957
    // The argument has a default, so the list always holds it; a value other than 0 stops the run
    // before any archive (`check_ncbi_options`).
    inputs.max_intron_length = args.max_intron_length.unwrap_or(0);
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(args.formatting.recorded_hitlist_start(), |size| size as i32);
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with TBLASTX's recorded options).
pub(crate) fn handle_for(
    args: &TblastxArgs,
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

/// The requests of a run's batches (`RunRequests`): nucleotide queries and subjects, the
/// `-subject` records `subjects` (`None` for `-db`), and with `-lcase_masking` (`lcase_masking`)
/// the lower-case runs of each query in both frames of a translated query.
pub(super) fn requests(subjects: Option<&[FastaRecord]>, lcase_masking: bool) -> RunRequests<'_> {
    RunRequests::new(
        RequestPolicy {
            query_protein: false,
            subject_protein: false,
            lcase: lcase_masking.then_some(LcaseFrames::PlusAndMinus),
        },
        subjects,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1718-1725
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
/// The `EffectiveSearchSpace` option of a batch whose first query has the statistics `first`
/// (`BatchArchive::searched_space`, `StrategyExport::note_searched`): the effective search space
/// of context 0 and that context's length in residues (`archive::RequestLength::space`), when NCBI
/// searched the batch.
pub(super) fn batch_searched_space(
    searched: bool,
    first: &TblastxQueryStats,
) -> Option<crate::report::archive::SearchedContext0> {
    searched.then_some(crate::report::archive::SearchedContext0 {
        space: first.context0_search_space,
        query_length: first.context0_aa_len,
    })
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
/// (`range` differs from `kTarget`, blast_aux.cpp:936-951), as BLASTX's.
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
                best: query.believed().map(|ids| ids.best.clone()),
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
                id: query.seq_id(),
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

/// Calls `each` with the queries and the archive of each batch of `data.spans`, in order
/// (`BatchArchive`: its request, `RunRequests::with_batch`, and its results).
fn for_each_batch(
    data: &ArchiveData<'_>,
    mut each: impl FnMut(Range<usize>, &BatchArchive<'_>) -> Result<()>,
) -> Result<()> {
    let requests = requests(Some(data.subject_records), data.lcase_masking);
    for span in data.spans.iter().cloned() {
        if span.is_empty() {
            continue;
        }
        let first = span.start;
        let searched = data.searched[first];
        let records = &data.queries[span.clone()];
        let placements = data.query_placements;

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
        requests.with_batch(records, placements.get(first), |request| {
            let batch = BatchArchive {
                request: *request,
                masks: &results_masks,
                ka_blocks,
                search_space,
                length_adjustment,
                searched_space: batch_searched_space(searched, stats),
                round: None,
            };
            each(span, &batch)
        })?;
    }
    Ok(())
}

/// Writes the archive of each batch of `data.spans` to `writer`; `warnings` writes the warnings
/// of a batch's queries (to the diagnostics, which saves them as the batch's messages) before
/// its archive, and its sink gets NCBI's posts of the bytes that the archive's stream replaces
/// (`write_archive`).
pub(super) fn write_tblastx_archives(
    handle: &ArchiveHandle,
    data: &ArchiveData<'_>,
    warnings: &mut QueryWarnings<'_>,
    writer: &mut dyn Write,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Tblastx,
        db_title: match data.subjects {
            SeqalignSubjects::Records { .. } => None,
            SeqalignSubjects::Database { title, .. } => Some(title),
        },
        gapped: false,
        program_type: None,
    };
    for_each_batch(data, |span, batch| {
        for index in span.clone() {
            warnings.before_query(index, &mut *writer)?;
        }
        let mut aligns = TblastxAligns {
            run: RunRecord {
                program: run.program,
                db_title: run.db_title,
                gapped: run.gapped,
                program_type: run.program_type,
            },
            data,
            span,
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

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:106-111
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3800-3837
/// tblastx's part of `-import_search_strategy` ([`import::apply`]): the arguments its command
/// line keeps and its table of the strategy's setters over its own options. tblastx searches
/// without gaps: the gap costs of its handle and the ungapped mode are its defaults, which the
/// table takes without a change.
impl import::StrategyArgs for TblastxArgs {
    const PROGRAM: Program = Program::Tblastx;

    fn runs(program: EProgram) -> bool {
        program == EProgram::Tblastx
    }

    fn imports(&self) -> bool {
        self.import_search_strategy.is_some() && self.import.is_none()
    }

    fn defaults() -> Result<Self> {
        match import::defaults_of("tblastx")? {
            crate::cli::Commands::Tblastx(args) => Ok(args),
            _ => unreachable!("the arguments of tblastx"),
        }
    }

    fn kept(self) -> import::Kept {
        let given = self.given.list().to_vec();
        import::Kept {
            path: self.import_search_strategy.expect("an import"),
            given,
            query: self.query,
            query_loc: self.query_loc,
            strand: Some(self.strand),
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
            mt_mode: None,
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
        if let Some(strand) = kept.strand {
            self.strand = strand;
        }
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

    fn set_task(&mut self, _task: &str) {}

    fn set(&mut self, setter: &Setter) -> Result<bool> {
        let program = Program::Tblastx.name();
        match setter {
            Setter::EvalueThreshold(v) => self.evalue = *v,
            Setter::WordSize(n) => self.word_size = Some(*n),
            Setter::WordThreshold(n) => self.threshold = Some(f64::from(*n)),
            Setter::WindowSize(n) => self.window_size = Some(*n),
            Setter::MatrixName(name) => self.matrix = Some(name.clone()),
            Setter::CullingLimit(n) => self.culling_limit = Some(*n),
            Setter::MaxHspsPerSubject(n) => {
                self.max_hsps = usize::try_from(*n).ok().filter(|n| *n > 0)
            }
            Setter::QueryCovHspPerc(v) => self.qcov_hsp_perc = Some(*v),
            Setter::MaskAtHash(on) => self.soft_masking = Some(*on),
            Setter::QueryGeneticCode(n) => {
                self.query_gencode = import::Start::genetic_code(*n, "query", program)?
            }
            Setter::DbGeneticCode(n) => {
                self.db_gencode = import::Start::genetic_code(*n, "database", program)?
            }
            Setter::SubjectBestHit => self.subject_besthit = true,
            Setter::BestHitOverhang(v) => self.best_hit_overhang = Some(*v),
            Setter::BestHitScoreEdge(v) => self.best_hit_score_edge = Some(*v),
            // A length other than 0 is NCBI's options error, as `-max_intron_length` gives it.
            Setter::LongestIntronLength(n) => self.max_intron_length = Some(*n),
            Setter::GapOpeningCost(_)
            | Setter::GapExtensionCost(_)
            | Setter::GappedMode(false)
            | Setter::SmithWatermanMode(false) => {}
            _ => return Ok(false),
        }
        Ok(true)
    }

    fn finish(&mut self, collected: import::Collected, _run: &import::ImportRun) -> Result<()> {
        let program = Program::Tblastx.name();
        collected.no_dust(program)?;
        collected.no_repeat(program)?;
        collected.no_nucleotide_algorithms(program)?;
        if let Some(value) = collected.seg.value() {
            self.seg = value;
        }
        // tblastx has no composition-based statistics: only mode 0 of a strategy runs.
        if collected
            .composition_value(false, program)?
            .is_some_and(|value| value != "0")
        {
            return Err(import::Start::unported(
                "a search strategy with composition-based statistics",
                program,
            ));
        }
        Ok(())
    }

    fn into_command(self) -> crate::cli::Commands {
        crate::cli::Commands::Tblastx(self)
    }
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod tests;
