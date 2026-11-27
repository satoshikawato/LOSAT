//! BLASTP's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:280-283
//! One archive per query batch, written after the batch's search; `report::archive` writes
//! the document, this module gathers what BLASTP holds for its other formats: the hit lists
//! (the Seq-aligns of outfmt 8), the SEG masks, and the Karlin-Altschul blocks and search
//! spaces of the batch's first query. BLASTP has one context per query, so the first valid
//! context of a batch is its first query (a first query without one gives both blocks the
//! `-1` triple), and the batch's `EffectiveSearchSpace` option is that query's search space.

use std::ops::Range;

use anyhow::Result;

use super::args::BlastpArgs;
use super::hsp::BlastpHitList;
use super::seqalign::{query_records, EngineCounts};
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::core::blast_util::ncbi_int4_from_double;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::OptionInputs;
use crate::report::archive::read::import;
use crate::report::archive::read::options::{EProgram, Setter};
use crate::report::archive::requests::{LcaseFrames, RequestPolicy, RunRequests};
use crate::report::archive::{
    ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive, KaBlock,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryNames, RunRecord, SubjectNames, SubjectRecord};
use crate::report::serial::seqalign::write_seq_aligns;
use crate::stats::tables::KarlinParams;

/// What a BLASTP run's arguments give its search strategy's requests
/// (`strategy::StrategyExport::for_args`): the archive's recorded options, the counts
/// of the `HitlistSize` and the database of `-db`.
impl crate::report::archive::strategy::StrategyArgs for BlastpArgs {
    fn strategy_inputs(&self) -> crate::report::archive::strategy::StrategyInputs {
        crate::report::archive::strategy::StrategyInputs {
            inputs: option_inputs(self),
            counts: crate::report::archive::strategy::HitlistArgs {
                max_target_seqs: self.max_target_seqs,
                num_descriptions: self.num_descriptions,
                num_alignments: self.num_alignments,
            },
            database: self.db.as_deref().map(|v| {
                crate::seqdb::db_argument_name_bytes(&crate::utils::os_bytes::os_str_bytes(v))
                    .to_vec()
            }),
        }
    }
}

/// The recorded options of a BLASTP run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2). An application whose search is BLASTP's (psiblast,
/// rpsblast) gives its own (`BlastpArgs::request_inputs`), made from these and its arguments.
pub(crate) fn option_inputs(args: &BlastpArgs) -> OptionInputs {
    // The archive's `HitlistSize` is that of outfmt 11, whoever made the other options.
    let hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(args.formatting.recorded_hitlist_start(), |size| size as i32);
    if let Some(inputs) = &args.request_inputs {
        let mut inputs = OptionInputs::clone(inputs);
        inputs.hitlist_size = hitlist_size;
        return inputs;
    }
    let mut inputs = OptionInputs::new(Program::Blastp, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    inputs.evalue = args.evalue;
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.matrix = args
        .matrix
        .as_deref()
        .map(|v| crate::utils::os_bytes::os_str_bytes(v).into_owned());
    // `-xdrop_gap` and `-xdrop_gap_final` as given (`SetGapXDropoff`, `SetGapXDropoffFinal`);
    // `-xdrop_ungap` has no entry (NCBI records none, `tblastn/archive.rs` `opt_xdropu`).
    inputs.xdrop_gap = args.xdrop_gap;
    inputs.xdrop_gap_final = args.xdrop_gap_final;
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps_per_subject.map(|count| count as i32);
    inputs.culling_limit = args.culling_limit;
    inputs.best_hit_overhang = args.best_hit_overhang;
    inputs.best_hit_score_edge = args.best_hit_score_edge;
    inputs.subject_besthit = args.subject_besthit;
    // `-threshold` is read as a double and recorded as an `int` (`static_cast<int>`), converted
    // as the x86-64 build does (`INT_MIN` past `Int4`, as the lookup table takes it:
    // `ncbi_int4_from_double`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:884-892
    inputs.threshold = args.threshold.map(ncbi_int4_from_double);
    inputs.window_size = args.window_size;
    inputs.ungapped = args.ungapped;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
    // `MaskAtHash` is recorded with the value given (the class always sets it: false when absent).
    inputs.soft_masking = args.soft_masking;
    inputs.seg = args.seg.as_deref().and_then(|v| {
        crate::report::archive::options::seg_filter(crate::utils::os_bytes::os_str_bytes(v))
    });
    inputs.imported = args.import.as_ref().map(|run| run.recorded.clone());
    inputs.application = args.import.as_ref().map(|run| run.app.ncbi_name());
    inputs.hitlist_size = hitlist_size;
    if let Ok((mode, unified_p)) = crate::blastinput::app::parse_comp_based_stats(
        crate::utils::os_bytes::os_str_bytes(&args.comp_based_stats),
        true,
        false,
    ) {
        inputs.comp_based_stats = mode as i32;
        inputs.unified_p = unified_p;
    }
    inputs.use_sw_tback = args.use_sw_tback;
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with BLASTP's recorded options).
pub(crate) fn handle_for(
    args: &BlastpArgs,
    archive_requested: bool,
    options: &crate::report::archive::OptionsHandle,
) -> Option<ArchiveHandle> {
    if !archive_requested {
        return None;
    }
    // NCBI keeps the name without the white space at its ends (`NStr::TruncateSpaces`,
    // blast_args.cpp:2472-2476), and the request's `database` is that name.
    let database = args.db.as_deref().map(|v| {
        crate::seqdb::db_argument_name_bytes(&crate::utils::os_bytes::os_str_bytes(v)).to_vec()
    });
    Some(ArchiveHandle::for_program(
        &option_inputs(args),
        database.as_deref(),
        options,
    ))
}

/// The requests of a run's batches (`RunRequests`): protein queries and subjects (the `-subject`
/// records `subjects`, `None` for `-db`), and with `-lcase_masking` (`lcase_masking`) one `notset`
/// mask of each query's lower-case runs (SEG is never one); `pssm` is the PSSM that every request
/// holds in place of the queries (`RunRequests::with_pssm_query`): the search strategy's of
/// psiblast's `-in_pssm` as read (`PsiSearch::request`), the archive's of a PSSM round
/// (`PsiArchive::pssm`).
pub(crate) fn requests(
    subjects: Option<&[FastaRecord]>,
    lcase_masking: bool,
    pssm: Option<std::sync::Arc<crate::report::serial::pssm::PssmWithParameters>>,
) -> RunRequests<'_> {
    RunRequests::new(
        RequestPolicy {
            query_protein: true,
            subject_protein: true,
            lcase: lcase_masking.then_some(LcaseFrames::NotSet),
        },
        subjects,
    )
    .with_pssm_query(pssm)
}

pub(crate) use crate::report::archive::error_lines;

/// The queries of a call's batches and what their archives read (`for_each_batch`).
pub(crate) struct CallBatches<'a> {
    /// The searched queries of the call (their records keep the whole letters of a `-query_loc`
    /// cut).
    pub queries: &'a [FastaRecord],
    pub query_placements: &'a Placements,
    /// `-subject`: every subject record, whole (the Bioseqs of the request).
    pub subject_records: Option<&'a [FastaRecord]>,
    /// The SEG masks and, with `-lcase_masking`, the lower-case runs of each query that the
    /// results carry, in record coordinates (`blast_engine::result_query_masks`).
    pub masks: &'a [Vec<(i32, i32)>],
    /// `-lcase_masking`.
    pub lcase_masking: bool,
    /// Whether the query's batch was not searched (every query of it invalid).
    pub unsearched: &'a [bool],
    /// The queries of each NCBI batch of the call, in order.
    pub spans: &'a [Range<usize>],
}

/// What psiblast's archive of one round sets in its request (`WriteArchive` with the round's
/// number, `blast_engine::report_write::psi_archive`).
pub(crate) struct PsiArchive {
    pub round: crate::report::archive::ArchiveRound,
    /// The PSSM that the request of a PSSM search holds in place of the query, with the round's
    /// K/A written back; `None` for a FASTA query's first round (its Bioseq).
    pub pssm: Option<std::sync::Arc<crate::report::serial::pssm::PssmWithParameters>>,
    /// psiblast `-phi_pattern`'s rounds after round 1 (item PS8b): the context 0 that PHI-BLAST's
    /// round 1 left in the options handle of the archive, whose `EffectiveSearchSpace` option it
    /// gives in place of the round's own (`RoundOutcome::phi_handle`).
    pub handle_space: Option<crate::report::archive::SearchedContext0>,
}

/// What the archive of a call needs: the call's batches (the batches of one block, or of the
/// whole search) and their results.
pub(crate) struct ArchiveCall<'a> {
    pub handle: &'a ArchiveHandle,
    /// psiblast's round; `None` for BLASTP and RPS-BLAST.
    pub psi: Option<&'a PsiArchive>,
    /// `!args.ungapped`: the archive's Seq-aligns are those of outfmt 8 (`write_blastp_seqalign`).
    pub gapped: bool,
    /// The app's program type (`BlastpApp::program_type`, `RunRecord::program_type`).
    pub program_type: Option<crate::report::record::ProgramType>,
    pub batches: CallBatches<'a>,
    pub query_lengths: &'a [usize],
    /// The subject names of the Seq-aligns.
    pub subjects: SubjectNames<'a>,
    pub subject_placements: &'a Placements,
    pub hit_lists: &'a [Option<BlastpHitList>],
    pub engine_counts: &'a EngineCounts,
}

/// The values of a batch's first query that `search-stats`, `ka-blocks` and the
/// `EffectiveSearchSpace` option take.
pub(crate) struct FirstContext {
    /// The ungapped block and, for a gapped search, the gapped one; `None` when the first query
    /// has no valid context (both blocks are then missing).
    pub karlin: Option<(KarlinParams, Option<KarlinParams>)>,
    /// The effective search space and length adjustment of that context (0 without one).
    pub search_space: i64,
    pub length_adjustment: i64,
    /// The `EffectiveSearchSpace` option: that of the first query's context, `-searchsp` in an
    /// invalid one too.
    pub context0_space: i64,
    /// That context's query length, 0 for an invalid one (`archive::RequestLength::space`).
    pub context0_length: i64,
}

/// The Seq-aligns of one batch of a call (the elements of outfmt 8, without the Seq-annot).
struct BlastpAligns<'a> {
    run: RunRecord<'a>,
    queries: QueryNames<'a>,
    subjects: SubjectNames<'a>,
    subject_placements: &'a Placements,
    hit_lists: &'a [Option<BlastpHitList>],
    engine_counts: &'a EngineCounts,
    /// The index of the batch's first query in the call.
    first: usize,
}

impl<'a> AlignSource for BlastpAligns<'a> {
    fn write<E: crate::report::serial::encoder::Encoder>(
        &mut self,
        out: &mut ArchiveOut<'_, E>,
    ) -> Result<()> {
        let mut segments = crate::report::serial::seqalign::Segments::default();
        let mut subject_records: Vec<SubjectRecord<'a>> = Vec::new();
        let mut hsp_records: Vec<HspRecord<'a>> = Vec::new();
        let hit_lists: &'a [Option<BlastpHitList>] = self.hit_lists;
        for (offset, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let q_idx = self.first + offset;
            query_records(
                hit_list,
                q_idx,
                &self.subjects,
                self.subject_placements,
                self.engine_counts,
                self.queries.placements.offset(q_idx) as i32,
                &mut subject_records,
                &mut hsp_records,
            );
            let query = self.queries.record(q_idx);
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

/// Calls `each` with the queries and the archive of each batch of `call`, in order
/// (`BatchArchive`: its request, `RunRequests::with_batch`, and its results).
/// `first_context(first_query)` gives the values of a batch's first query; `psi` is psiblast's
/// round.
pub(crate) fn for_each_batch(
    call: &CallBatches<'_>,
    psi: Option<&PsiArchive>,
    first_context: &dyn Fn(usize) -> FirstContext,
    mut each: impl FnMut(Range<usize>, &BatchArchive<'_>) -> Result<()>,
) -> Result<()> {
    let whole = [0..call.queries.len()];
    let spans: &[Range<usize>] = if call.spans.is_empty() {
        &whole
    } else {
        call.spans
    };
    // The requests' queries are the call's records, or the PSSM of a PSSM search.
    let requests = requests(
        call.subject_records,
        call.lcase_masking,
        psi.and_then(|psi| psi.pssm.clone()),
    );
    for span in spans.iter().cloned() {
        if span.is_empty() {
            continue;
        }
        let first = span.start;
        let searched = !call.unsearched[first];
        let context = first_context(first);

        // The results' masks are those of the search, so only a searched batch has them: one
        // `notset` mask per query with the SEG intervals and lower-case runs (a protein query).
        // `call.masks` are the results' (`blast_engine::result_query_masks`: none for psiblast).
        let mut masks: Vec<MaskEntry> = Vec::new();
        if searched {
            for index in span.clone() {
                if !call.masks[index].is_empty() {
                    masks.push(MaskEntry {
                        id: call.queries[index].local_id.clone(),
                        best: call.queries[index].believed().map(|ids| ids.best.clone()),
                        frame: Frame::NotSet,
                        intervals: call.masks[index].clone(),
                    });
                }
            }
        }

        let ka_blocks = match &context.karlin {
            Some((ungapped, gapped)) => [
                ka_block(ungapped, false),
                gapped
                    .as_ref()
                    .map_or(KaBlock::missing(true), |gapped| ka_block(gapped, true)),
            ],
            None => [KaBlock::missing(false), KaBlock::missing(true)],
        };
        requests.with_batch(
            &call.queries[span.clone()],
            call.query_placements.get(first),
            |request| {
                let batch = BatchArchive {
                    request: *request,
                    masks: &masks,
                    ka_blocks,
                    search_space: context.search_space,
                    length_adjustment: context.length_adjustment,
                    searched_space: match psi.and_then(|psi| psi.handle_space) {
                        Some(handle) => Some(handle),
                        None => searched.then_some(crate::report::archive::SearchedContext0 {
                            space: context.context0_space,
                            query_length: context.context0_length,
                        }),
                    },
                    round: psi.map(|psi| psi.round),
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
pub(crate) fn write_blastp_archives(
    call: &ArchiveCall<'_>,
    warnings: &mut QueryWarnings<'_>,
    first_context: &dyn Fn(usize) -> FirstContext,
    writer: &mut dyn crate::report::ostream::ReportWrite,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Blastp,
        db_title: call.subjects.db_title(),
        // An ungapped search's Seq-aligns are Dense-diags, as in outfmt 8 (`seqalign.rs`).
        gapped: call.gapped,
        program_type: call.program_type,
    };
    for_each_batch(&call.batches, call.psi, first_context, |span, batch| {
        for index in span.clone() {
            warnings.before_query(index, &mut *writer)?;
        }
        let mut aligns = BlastpAligns {
            run,
            queries: QueryNames {
                records: call.batches.queries,
                searched_lengths: call.query_lengths,
                placements: call.batches.query_placements,
            },
            subjects: call.subjects,
            subject_placements: call.subject_placements,
            hit_lists: &call.hit_lists[span.clone()],
            engine_counts: call.engine_counts,
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

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:128-133
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3800-3837
/// blastp's part of `-import_search_strategy` ([`import::apply`]): the arguments its command
/// line keeps and its table of the strategy's setters over its own options.
impl import::StrategyArgs for BlastpArgs {
    const PROGRAM: Program = Program::Blastp;

    fn runs(program: EProgram) -> bool {
        program == EProgram::Blastp
    }

    fn imports(&self) -> bool {
        self.import_search_strategy.is_some()
            && self.import.is_none()
            && self.app == super::args::BlastpApp::Blastp
    }

    fn defaults() -> Result<Self> {
        match import::defaults_of("blastp")? {
            crate::cli::Commands::Blastp(args) => Ok(args),
            _ => unreachable!("the arguments of blastp"),
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
            Setter::GapTrigger(v) => self.gap_trigger = Some(*v),
            Setter::CullingLimit(n) => self.culling_limit = Some(*n),
            Setter::MaxHspsPerSubject(n) => {
                self.max_hsps_per_subject = usize::try_from(*n).ok().filter(|n| *n > 0)
            }
            Setter::QueryCovHspPerc(v) => self.qcov_hsp_perc = Some(*v),
            Setter::SubjectBestHit => self.subject_besthit = true,
            Setter::BestHitOverhang(v) => self.best_hit_overhang = Some(*v),
            Setter::BestHitScoreEdge(v) => self.best_hit_score_edge = Some(*v),
            Setter::MaskAtHash(on) => self.soft_masking = Some(*on),
            Setter::SmithWatermanMode(on) => self.use_sw_tback = *on,
            Setter::GappedMode(gapped) => self.ungapped = !*gapped,
            Setter::SumStatisticsMode(on) => self.sum_stats = *on,
            _ => return Ok(false),
        }
        Ok(true)
    }

    fn finish(&mut self, collected: import::Collected, _run: &import::ImportRun) -> Result<()> {
        let program = Program::Blastp.name();
        collected.no_dust(program)?;
        collected.no_repeat(program)?;
        collected.no_nucleotide_algorithms(program)?;
        if let Some(value) = collected.seg.value() {
            self.seg = Some(value.into());
        }
        if let Some(value) = collected.composition_value(true, program)? {
            self.comp_based_stats = value.into();
        }
        Ok(())
    }

    fn into_command(self) -> crate::cli::Commands {
        crate::cli::Commands::Blastp(self)
    }
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod tests;
