//! RPSTBLASTN's search of BLASTX's query batches (series PR, item PS12a, DECISIONS 38 D1): the six
//! frames of each nucleotide query, set up as BLASTX sets them up (`prepare_unit_queries`: frames,
//! strands, genetic codes, masks, `-query_loc`; `score_block` and `effective_lengths`: the frames'
//! Karlin-Altschul blocks and search spaces), are the "subjects" of the reversed RPS search; its
//! final HSPs become BLASTX's results (`blastx_lists`), which BLASTX's reports write. NCBI's
//! `eBlastTypeRpsTblastn` is `PSSM_SUBJECT_MASK | eBlastTypeBlastx`.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_program.h:85
//! ```c
//!     eBlastTypeRpsTblastn    = (PSSM_SUBJECT_MASK | eBlastTypeBlastx),
//! ```
//! NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:191-209
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//!
//!             SaveSearchStrategy(args, m_CmdLineArgs, queries, opts_hndl);
//!
//!             CRef<CSearchResultSet> results;
//!
//!             if (m_CmdLineArgs->ExecuteRemotely()) {
//!                 ...
//!             } else {
//!             	CLocalRPSBlast  local_search (query_batch, db_args->GetDatabaseName(), opts_hndl, args[kArgNumThreads].AsInteger() );
//!             	results = local_search.Run();
//!             }
//! ```

use std::io::Write;

use anyhow::{bail, Result};

use super::engine::RpsSearch;
use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHsp};
use crate::algorithm::blastx::args::{ResolvedOptions, RpsTblastnSearch};
use crate::algorithm::blastx::database::BlastxDatabase;
use crate::algorithm::blastx::input::FastaRecord;
use crate::algorithm::blastx::native::{engine, BatchRender, QueryBatchReader};
use crate::algorithm::blastx::parameters::{effective_lengths, score_block, ContextParameters};
use crate::algorithm::blastx::preliminary::PreliminaryHsp;
use crate::algorithm::blastx::query_setup::{prepare_unit_queries, PreparedQueryBatch};
use crate::algorithm::blastx::results::{Hsp, HspList};
use crate::algorithm::blastx::runtime::BatchResults;
use crate::blastinput::report_counts::FormattingArgs;
use crate::seqdb::rps::{RpsOpenError, RpsVolumeFiles};
use crate::stats::KarlinParams;

/// The six contexts of a translated query (`NUM_FRAMES`).
const NUM_FRAMES: usize = 6;

/// The RPS search of an RPSTBLASTN run: `RpsSearch` (rpsblast's options of the RPS engine), the
/// formatting options that BLASTX does not have and `-ungapped`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TranslatedRpsSearch {
    pub(crate) search: RpsSearch,
    pub(crate) formatting: FormattingArgs,
    /// `-ungapped` (with composition-based statistics 0: NCBI's options error at the others is
    /// BLASTX's `resolve`).
    pub(crate) ungapped: bool,
}

impl RpsTblastnSearch for TranslatedRpsSearch {
    fn run_query_batches(
        &self,
        options: &ResolvedOptions,
        database: Option<&BlastxDatabase>,
        read: &mut QueryBatchReader<'_>,
        diagnostics: &mut dyn Write,
        search_error: &dyn Fn(anyhow::Error) -> anyhow::Error,
        render: &mut BatchRender<'_>,
    ) -> Result<usize> {
        let Some(database) = database else {
            bail!("a search without -db is not supported by LOSAT's RPSTBLASTN");
        };
        let mut processed = 0usize;
        let diagnostics = std::cell::RefCell::new(diagnostics);
        let volumes = SortedVolumes::of(database);
        let (db_length, db_count) = database_lengths(database);
        read(
            &mut |queries| {
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
                // ```c++
                // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
                //     : m_QueryVector(& queries)
                // {
                //     if (queries.Empty()) {
                //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
                //     }
                // ```
                if queries.is_empty() {
                    return Err(engine("Empty CBlastQueryVector"));
                }
                let mut results = vec![self
                    .search_batch(queries, options, &volumes, db_length, db_count)
                    .map_err(|error| match error.downcast::<RpsBatchError>() {
                        // NCBI's errors of the set-up of the batch (`BLAST engine error`).
                        Ok(RpsBatchError::Setup(error)) => search_error(error),
                        Ok(RpsBatchError::Engine(text)) => engine(text),
                        Err(refusal) => refusal,
                    })?];
                render(
                    &mut **diagnostics.borrow_mut(),
                    queries,
                    &mut results,
                    processed,
                )?;
                processed += queries.len();
                Ok(())
            },
            &mut |bytes: &[u8]| Ok(diagnostics.borrow_mut().write_all(bytes)?),
        )?;
        Ok(processed)
    }

    fn formatting(&self) -> FormattingArgs {
        self.formatting
    }
}

/// What ends the search of a query batch besides LOSAT's refusals: NCBI's error of the batch's
/// set-up (BLASTX's, `prepare_unit_queries`), or NCBI's engine error with its text.
#[derive(Debug)]
enum RpsBatchError {
    Setup(anyhow::Error),
    Engine(String),
}

impl std::fmt::Display for RpsBatchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Setup(error) => write!(f, "{error}"),
            Self::Engine(text) => f.write_str(text),
        }
    }
}

impl std::error::Error for RpsBatchError {}

/// The volumes of the profile database in NCBI's search order (`SeqDB_CompareVolume`: by base
/// name, then path; rpsblast_local.cpp's sorted volume list), with the first ordinal id of each.
struct SortedVolumes<'a> {
    database: &'a BlastxDatabase,
    /// The index of each volume in the database's order, sorted.
    order: Vec<usize>,
    /// The first ordinal id of each volume, in the database's order.
    starts: Vec<u32>,
}

impl<'a> SortedVolumes<'a> {
    fn of(database: &'a BlastxDatabase) -> Self {
        let volumes = database.x_seqdb().volumes();
        let mut starts = Vec::with_capacity(volumes.len());
        let mut next = 0u32;
        for volume in volumes {
            starts.push(next);
            next += volume.num_oids();
        }
        let mut order: Vec<usize> = (0..volumes.len()).collect();
        order.sort_by(|&a, &b| {
            crate::seqdb::SeqDb::compare_volume_paths(volumes[a].base(), volumes[b].base())
        });
        Self {
            database,
            order,
            starts,
        }
    }

    /// The RPS files of the volume `index` (database order), mapped: NCBI's errors of a missing or
    /// malformed `.rps`, `.loo` or `.aux` (rc 3, `RpsBatchError::Engine`), LOSAT's refusals of the
    /// files and of a matrix or gap costs other than BLOSUM62 11/1 (the `.aux`'s), whose
    /// statistics LOSAT's tables do not hold for the RPS search.
    fn open(&self, index: usize) -> Result<RpsVolumeFiles> {
        let volume = &self.database.x_seqdb().volumes()[index];
        let lengths = (0..volume.num_oids())
            .map(|oid| volume.protein_length(oid))
            .collect::<std::result::Result<Vec<u32>, _>>()
            .map_err(|error| {
                crate::algorithm::common::database::database_error(error, "RPSTBLASTN")
            })?;
        let files = RpsVolumeFiles::open(volume.base(), &lengths).map_err(|error| match error {
            RpsOpenError::Ncbi(text) => anyhow::Error::new(RpsBatchError::Engine(text)),
            RpsOpenError::Refused(text) => anyhow::anyhow!(text),
        })?;
        let aux = &files.aux;
        if aux.matrix != "BLOSUM62" || aux.gap_open != 11 || aux.gap_extend != 1 {
            bail!(
                "a profile database of the matrix {} with gap costs {}/{} is not supported by LOSAT's RPSTBLASTN",
                aux.matrix,
                aux.gap_open,
                aux.gap_extend
            );
        }
        Ok(files)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-720
// ```c
//    /* use overriding value from effective lengths options or the real value
//       from effective lengths parameters. */
//    if (eff_len_options->db_length > 0)
//       db_length = eff_len_options->db_length;
//    else
// ```
/// The database's length and number of sequences for the effective lengths of the frames: the
/// whole profile database's (every volume's search reads the totals, `x_AdjustDbSize` without
/// `-dbsize`); the database is protein, so the length is not divided by 3.
fn database_lengths(database: &BlastxDatabase) -> (i64, i64) {
    let pairwise = database.pairwise_database();
    (pairwise.total_letters as i64, pairwise.num_sequences as i64)
}

impl TranslatedRpsSearch {
    /// The search of one query batch: BLASTX's set-up of its frames (`prepare_unit_queries`,
    /// `score_block`, `effective_lengths`), the RPS search of every volume and the batch's
    /// results. A batch without a valid frame is not searched (NCBI's `CheckInternalData`): its
    /// queries get the Karlin-Altschul warning of BLASTX's reports (`search_skipped`).
    fn search_batch(
        &self,
        queries: &[FastaRecord],
        options: &ResolvedOptions,
        volumes: &SortedVolumes<'_>,
        db_length: i64,
        db_count: i64,
    ) -> Result<BatchResults> {
        let prepared =
            prepare_unit_queries(queries, None, options).map_err(RpsBatchError::Setup)?;
        // BLOSUM62 (the `.aux`'s matrix, `SortedVolumes::open`): BLASTX's ideal block, which
        // `check_ideal` takes for RPSTBLASTN as for BLASTX.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2746-2748
        // ```c
        //     Boolean check_ideal =
        //        (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
        //         program == eBlastTypeRpsTblastn);
        // ```
        let mut parameters = score_block(&prepared);
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181
        // ```c++
        //
        //     int status = m_PrelimSearch->CheckInternalData();
        //     if (status != 0)
        //     {
        //          // Search was not run, but we send back an empty CSearchResultSet.
        //          CRef<ILocalQueryData> local_query_data = m_QueryFactory->MakeLocalQueryData(m_Opts);
        // ```
        if !parameters.iter().any(|p| p.valid) {
            return Ok(unsearched_batch(prepared, parameters));
        }
        effective_lengths(
            &prepared,
            &mut parameters,
            options,
            db_length,
            db_count,
            None,
        )
        .map_err(RpsBatchError::Setup)?;
        let frames = frame_queries(&prepared, &parameters);
        // Volume after volume, each mapped for the batch's search alone, as NCBI's
        // `CLocalRPSBlast` searches each query batch.
        let mut hit_lists: Vec<Option<BlastpHitList>> = vec![None; frames.len()];
        for (position, &index) in volumes.order.iter().enumerate() {
            let files = volumes.open(index)?;
            if position == 0 {
                self.check_search()?;
            } else {
                bail!("a profile database of several volumes is not supported by LOSAT's RPSTBLASTN yet");
            }
            let mut lists = search_volume(&self.search, &files, &frames)?;
            drop(files);
            // The volume's lists named by the database's ordinal ids (`rps_hooks`'
            // `shift_ordinal_ids`).
            let start = volumes.starts[index];
            for list in lists.iter_mut().flatten() {
                let count = list.hsplist_count;
                for hsp_list in list.hsplist_array.iter_mut().take(count) {
                    hsp_list.oid += start;
                    for hsp in &mut hsp_list.hsps {
                        hsp.s_idx += start;
                    }
                }
            }
            hit_lists = lists;
        }
        Ok(BatchResults {
            query_ordinal: 0,
            queries: hit_lists
                .into_iter()
                .enumerate()
                .map(|(query, list)| blastx_lists(list, query))
                .collect(),
            prepared,
            parameters,
            search_skipped: false,
            gapped_k: None,
        })
    }

    /// LOSAT's refusals where NCBI's search starts, after NCBI's errors of the options, the
    /// queries and the database's files: what NCBI BLAST+ 2.17.0's rpstblastn does there has no
    /// result that can be checked.
    fn check_search(&self) -> Result<()> {
        if self.search.sum_stats {
            bail!(
                "-sum_stats true (NCBI BLAST+ 2.17.0's rpstblastn links HSPs with the profiles' indexes as query indexes and reads past its arrays: its output depends on memory contents or it crashes) is not supported by LOSAT's RPSTBLASTN"
            );
        }
        // DECISIONS 38 D3: every probed query (`out/u.*` of the design) segfaults after the
        // search, valgrind: a read at address 0 in `BLAST_ComputeTraceback_MT`.
        if self.ungapped {
            bail!(
                "-ungapped (NCBI BLAST+ 2.17.0's rpstblastn crashes after an ungapped search) is not supported by LOSAT's RPSTBLASTN"
            );
        }
        Ok(())
    }
}

/// The results of a batch that is not searched: no list for any query, and BLASTX's reports
/// warn for each (`report::query_warning`).
fn unsearched_batch(
    prepared: PreparedQueryBatch,
    parameters: Vec<ContextParameters>,
) -> BatchResults {
    BatchResults {
        query_ordinal: 0,
        queries: vec![Vec::new(); prepared.original_lengths.len()],
        prepared,
        parameters,
        search_skipped: true,
        gapped_k: None,
    }
}

/// The RPS search of one volume for the frames of a batch's queries: each query's hit list, in
/// the batch's order. The frames of a query are the subjects of the reversed search (NCBI's
/// frame loop of `s_RPSPreliminarySearchEngine`, blast_engine.c:757-846); the frames' results
/// combine by score before the traceback, which aligns each HSP on its frame.
fn search_volume(
    _search: &RpsSearch,
    _files: &RpsVolumeFiles,
    _frames: &[RpsFrameQuery<'_>],
) -> Result<Vec<Option<BlastpHitList>>> {
    bail!("the RPS search of translated queries is not supported by LOSAT's RPSTBLASTN yet")
}

/// One frame (context) of a user query as the RPS search reads it: NCBI's `subject->frame`,
/// `subject->sequence` and `subject->length` of the frame loop, and the frame's statistics.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-813
/// ```c
///     for (context=first_context; context<=last_context; context++) {
///         BlastHSPList* hsp_list_for_chunks = NULL;
///
///         if (kTranslatedSubject) {
///             Uint4 i;
///             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
///             subject->sequence = translation_buffer + frame_offsets[context] + 1;
///             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
///             if (subject->length > 0) stat_length = subject->length;
/// ```
#[derive(Debug, Clone, Copy)]
#[allow(
    dead_code,
    reason = "the frames' search (`search_volume`) reads them once the engine runs frames (item PS12a's second part)"
)]
pub(crate) struct RpsFrame<'a> {
    /// +1, +2, +3, -1, -2 or -3.
    pub(crate) frame: i8,
    /// The frame's letters between two sentinels as searched (SEG and lower-case masks); empty
    /// for a frame without letters (the other strand of `-strand`, a query shorter than the
    /// frame), which the search does not scan.
    pub(crate) aa_seq: &'a [u8],
    /// The frame's letters as read (between two sentinels) where the searched ones differ: the
    /// reports' identities count on them, as rpsblast's (`traceback::formatter_identities`).
    pub(crate) aa_seq_nomask: Option<&'a [u8]>,
    /// Whether the frame's Karlin-Altschul block exists (an invalid frame gets no HSP).
    pub(crate) valid: bool,
    /// `kbp_std[context]`: the frame's composition with BLOSUM62, the ideal block where its
    /// lambda is the larger (`check_ideal`).
    pub(crate) kbp_std: KarlinParams,
    /// `eff_searchsp` of the context.
    pub(crate) search_space: i64,
}

#[allow(
    dead_code,
    reason = "the frames' search (`search_volume`) reads them once the engine runs frames (item PS12a's second part)"
)]
impl RpsFrame<'_> {
    /// The frame's letters without the sentinels (`subject->length`).
    pub(crate) fn length(&self) -> usize {
        self.aa_seq.len().saturating_sub(2)
    }
}

/// The six frames of one user query, in NCBI's context order (+1, +2, +3, -1, -2, -3).
#[derive(Debug, Clone)]
#[allow(
    dead_code,
    reason = "the frames' search (`search_volume`) reads them once the engine runs frames (item PS12a's second part)"
)]
pub(crate) struct RpsFrameQuery<'a> {
    pub(crate) frames: [RpsFrame<'a>; NUM_FRAMES],
}

#[allow(
    dead_code,
    reason = "the frames' search (`search_volume`) reads them once the engine runs frames (item PS12a's second part)"
)]
impl RpsFrameQuery<'_> {
    /// `stat_length` of the preliminary e-values: the length of the last frame that has letters
    /// (frame -3 with both strands, +3 with `-strand plus`), not the query's length; 0 when no
    /// frame has letters.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:813
    /// ```c
    ///             if (subject->length > 0) stat_length = subject->length;
    /// ```
    pub(crate) fn prelim_stat_length(&self) -> usize {
        self.frames
            .iter()
            .rev()
            .map(RpsFrame::length)
            .find(|&length| length > 0)
            .unwrap_or(0)
    }
}

/// The frames of each query of a prepared batch with the statistics of its contexts
/// (`score_block`, `effective_lengths`): BLASTX's layout of `sequence_start` holds each frame
/// with letters between two sentinels at its context's offset.
pub(crate) fn frame_queries<'a>(
    prepared: &'a PreparedQueryBatch,
    parameters: &[ContextParameters],
) -> Vec<RpsFrameQuery<'a>> {
    let frame = |context: usize| -> RpsFrame<'a> {
        let c = &prepared.contexts[context];
        let p = &parameters[context];
        let letters = |buffer: &'a [u8]| -> &'a [u8] {
            if c.length == 0 {
                &[]
            } else {
                &buffer[c.offset..c.offset + c.length + 2]
            }
        };
        let aa_seq = letters(&prepared.sequence_start);
        let nomask = letters(&prepared.sequence_start_nomask);
        RpsFrame {
            frame: c.frame,
            aa_seq,
            aa_seq_nomask: (nomask != aa_seq).then_some(nomask),
            valid: p.valid && c.length > 0,
            kbp_std: p.ungapped,
            search_space: p.search_space,
        }
    };
    (0..prepared.original_lengths.len())
        .map(|query| RpsFrameQuery {
            frames: std::array::from_fn(|i| frame(query * NUM_FRAMES + i)),
        })
        .collect()
}

/// BLASTX's context of a frame within its query (`BLAST_FrameToContext` of blastx): +1, +2, +3
/// are 0, 1, 2; -1, -2, -3 are 3, 4, 5.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1211-1220
/// ```c
/// Int4 BLAST_FrameToContext(Int2 frame, EBlastProgramType program)
/// {
///     if (Blast_QueryIsTranslated(program) ||
///         Blast_SubjectIsTranslated(program)) {
///         ASSERT(frame >= -3 && frame <= 3 && frame != 0);
///         if (frame > 0) {
///             return frame - 1;
///         } else {
///             return 2 - frame;
///         }
/// ```
pub(crate) fn frame_context(frame: i32) -> usize {
    debug_assert!((-3..=3).contains(&frame) && frame != 0);
    if frame > 0 {
        (frame - 1) as usize
    } else {
        (2 - frame) as usize
    }
}

/// An RPS final HSP as BLASTX's results hold it: the roles are back (the query side is the
/// user query's frame, in protein letters; the subject side the profile's consensus), so the
/// coordinates become BLASTX's 0-based half-open ones on the frame and the profile, the context
/// is the query's frame context, and the scores, the e-value, the counts, the edit script and
/// the composition-based method are the traceback's.
fn blastx_hsp(hsp: &BlastpHsp, query_index: usize) -> Hsp {
    let to_i32 = |value: usize| i32::try_from(value).expect("a position below 2^31");
    Hsp {
        hsp: PreliminaryHsp {
            context: query_index * NUM_FRAMES + frame_context(hsp.query_frame),
            frame: hsp.query_frame as i8,
            score: hsp.raw_score,
            q_start: to_i32(hsp.q_start) - 1,
            q_end: to_i32(hsp.q_end),
            q_gapped_start: hsp.gapped_q_start,
            s_start: to_i32(hsp.s_start) - 1,
            s_end: to_i32(hsp.s_end),
            s_gapped_start: hsp.gapped_s_start,
        },
        // No sum statistics (`do_sum_stats` is FALSE for RPSTBLASTN): no linked set.
        num: 0,
        evalue: hsp.e_value,
        bit_score: hsp.bit_score,
        identity: hsp.num_ident,
        positive: hsp.num_positives,
        edit_script: hsp.gap_info.clone().unwrap_or_default(),
        composition_method: i32::from(hsp.comp_adjustment_method),
    }
}

/// A query's final hit list as BLASTX's results: one list per profile in the hit list's order,
/// each with its HSPs in their order, named by the profile's ordinal id in the database.
fn blastx_lists(hit_list: Option<BlastpHitList>, query_index: usize) -> Vec<HspList> {
    let Some(hit_list) = hit_list else {
        return Vec::new();
    };
    hit_list
        .hsplist_array
        .into_iter()
        .take(hit_list.hsplist_count)
        .filter(|list| !list.hsps.is_empty())
        .map(|list| {
            let mut converted = HspList {
                query_index,
                oid: list.oid as usize,
                hsps: list
                    .hsps
                    .iter()
                    .map(|hsp| blastx_hsp(hsp, query_index))
                    .collect(),
                best_evalue: list.best_evalue,
            };
            converted.refresh();
            converted
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::hsp::BlastpHspList;
    use crate::common::GapEditOp;

    fn options(extra: &[&str]) -> ResolvedOptions {
        let mut argv = vec!["losat", "rpstblastn", "-query", "q.fna", "-db", "cdd"];
        argv.extend_from_slice(&["-comp_based_stats", "0"]);
        argv.extend_from_slice(extra);
        let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).unwrap();
        let crate::cli::Commands::Rpstblastn(args) = cli.command else {
            unreachable!()
        };
        args.into_blastx_args().unwrap().resolve().unwrap()
    }

    fn record(id: &str, letters: &[u8]) -> FastaRecord {
        FastaRecord {
            internal_id: id.to_string(),
            title: String::new(),
            sequence: letters.to_vec(),
            lowercase_masks: Vec::new(),
            warnings: Vec::new(),
            source: None,
        }
    }

    /// 60 letters of a coding sequence (20 codons) whose six frames all have letters.
    const DNA: &[u8] = b"ATGGCTAAAGGTGAACGTCTGATCGCTGTTGCTAAAGAACTGGGTCTGACCGCTAAAGAA";

    // blast_engine.c:805-813 (quoted at `RpsFrame`): the frames of each query in NCBI's context
    // order, their letters at BLASTX's offsets between two sentinels, and the frames of the strand
    // that `-strand` leaves out without letters.
    #[test]
    fn the_frames_of_a_query_are_blastxs_contexts() {
        let options = options(&[]);
        let records = [record("a", DNA), record("b", &DNA[..31])];
        let prepared = prepare_unit_queries(&records, None, &options).unwrap();
        let parameters = score_block(&prepared);
        let queries = frame_queries(&prepared, &parameters);
        assert_eq!(queries.len(), 2);
        let frames: Vec<i8> = queries[0].frames.iter().map(|f| f.frame).collect();
        assert_eq!(frames, [1, 2, 3, -1, -2, -3]);
        let lengths: Vec<usize> = queries[0].frames.iter().map(RpsFrame::length).collect();
        assert_eq!(lengths, [20, 19, 19, 20, 19, 19]);
        let lengths: Vec<usize> = queries[1].frames.iter().map(RpsFrame::length).collect();
        assert_eq!(lengths, [10, 10, 9, 10, 10, 9]);
        for (query, frame_query) in queries.iter().enumerate() {
            for (i, frame) in frame_query.frames.iter().enumerate() {
                let context = &prepared.contexts[query * NUM_FRAMES + i];
                assert_eq!(frame.aa_seq.first(), Some(&0));
                assert_eq!(frame.aa_seq.last(), Some(&0));
                assert_eq!(
                    &frame.aa_seq[1..frame.aa_seq.len() - 1],
                    &prepared.sequence_start
                        [context.offset + 1..context.offset + 1 + context.length]
                );
                assert!(frame.valid);
                assert_eq!(frame.aa_seq_nomask, None);
                assert_eq!(frame_context(i32::from(frame.frame)), i);
            }
        }
        // The first frame translates the coding sequence (BLASTX's genetic code 1).
        let translated: Vec<u8> = queries[0].frames[0].aa_seq[1..21]
            .iter()
            .map(|&letter| b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ"[usize::from(letter)])
            .collect();
        assert_eq!(translated, b"MAKGERLIAVAKELGLTAKE");
        let plus = options_and_frames(&["-strand", "plus"]);
        assert_eq!(plus, [20, 19, 19, 0, 0, 0]);
        let minus = options_and_frames(&["-strand", "minus"]);
        assert_eq!(minus, [0, 0, 0, 20, 19, 19]);
    }

    fn options_and_frames(extra: &[&str]) -> Vec<usize> {
        let options = options(extra);
        let records = [record("a", DNA)];
        let prepared = prepare_unit_queries(&records, None, &options).unwrap();
        let parameters = score_block(&prepared);
        let queries = frame_queries(&prepared, &parameters);
        for frame in &queries[0].frames {
            assert_eq!(frame.valid, frame.length() > 0);
        }
        queries[0].frames.iter().map(RpsFrame::length).collect()
    }

    // blast_engine.c:813 (quoted at `prelim_stat_length`): the last frame with letters.
    #[test]
    fn the_preliminary_stat_length_is_the_last_frame_with_letters() {
        let options_both = options(&[]);
        let records = [record("a", &DNA[..59])];
        let prepared = prepare_unit_queries(&records, None, &options_both).unwrap();
        let parameters = score_block(&prepared);
        let queries = frame_queries(&prepared, &parameters);
        // Frames of 59 letters: +1 19, +2 19, +3 19, -1 19, -2 19, -3 19 ((59 - i) / 3).
        assert_eq!(
            queries[0].prelim_stat_length(),
            queries[0].frames[5].length()
        );
        let options_plus = options(&["-strand", "plus"]);
        let prepared = prepare_unit_queries(&[record("b", DNA)], None, &options_plus).unwrap();
        let parameters = score_block(&prepared);
        let queries = frame_queries(&prepared, &parameters);
        assert_eq!(queries[0].prelim_stat_length(), 19);
        // A query of 4 letters: one letter in frames +1, -1 only (frames +2, +3: (4 - 1) / 3 =
        // 1, (4 - 2) / 3 = 0 letters).
        let prepared = prepare_unit_queries(&[record("c", b"ACGT")], None, &options_both).unwrap();
        let parameters = score_block(&prepared);
        let queries = frame_queries(&prepared, &parameters);
        let lengths: Vec<usize> = queries[0].frames.iter().map(RpsFrame::length).collect();
        assert_eq!(lengths, [1, 1, 0, 1, 1, 0]);
        assert_eq!(queries[0].prelim_stat_length(), 1);
    }

    // blast_util.c:1211-1220 (quoted at `frame_context`).
    #[test]
    fn frames_have_blastxs_contexts() {
        let contexts: Vec<usize> = [1, 2, 3, -1, -2, -3].map(frame_context).to_vec();
        assert_eq!(contexts, [0, 1, 2, 3, 4, 5]);
    }

    fn hsp(frame: i32, score: i32, q: (usize, usize), s: (usize, usize)) -> BlastpHsp {
        BlastpHsp {
            identity: 50.0,
            length: q.1 - q.0 + 2,
            mismatch: 1,
            gapopen: 1,
            q_start: q.0,
            q_end: q.1,
            s_start: s.0,
            s_end: s.1,
            gapped_q_start: (q.0 + 2) as i32,
            gapped_s_start: (s.0 + 2) as i32,
            e_value: 1e-5,
            bit_score: 30.5,
            num_ident: 7,
            query_context: 0,
            query_frame: frame,
            subject_frame: 0,
            query_length: 100,
            q_idx: 0,
            s_idx: 0,
            raw_score: score,
            // The query's columns (1-based closed), with one deletion after the fifth.
            gap_info: Some(match (q.1 - q.0 + 1) as u32 {
                columns if columns > 5 => vec![
                    GapEditOp::Sub(5),
                    GapEditOp::Del(1),
                    GapEditOp::Sub(columns - 5),
                ],
                columns => vec![GapEditOp::Sub(columns)],
            }),
            num_positives: 9,
            comp_adjustment_method: 0,
        }
    }

    // The conversion of the RPS traceback's final HSPs (1-based closed coordinates after the role
    // swap, `traceback::swap_roles`) into BLASTX's (0-based half-open, context = query * 6 +
    // the frame's context).
    #[test]
    fn rps_hsps_become_blastx_hsps() {
        let converted = blastx_hsp(&hsp(-2, 41, (11, 30), (3, 23)), 2);
        assert_eq!(
            converted.hsp,
            PreliminaryHsp {
                context: 2 * 6 + 4,
                frame: -2,
                score: 41,
                q_start: 10,
                q_end: 30,
                q_gapped_start: 13,
                s_start: 2,
                s_end: 23,
                s_gapped_start: 5,
            }
        );
        assert_eq!(
            (
                converted.evalue,
                converted.bit_score,
                converted.identity,
                converted.positive,
                converted.num,
                converted.composition_method
            ),
            (1e-5, 30.5, 7, 9, 0, 0)
        );
        assert_eq!(
            converted.edit_script,
            [GapEditOp::Sub(5), GapEditOp::Del(1), GapEditOp::Sub(15)]
        );
        // Every frame's context.
        for (frame, context) in [(1, 0), (2, 1), (3, 2), (-1, 3), (-2, 4), (-3, 5)] {
            assert_eq!(
                blastx_hsp(&hsp(frame, 1, (1, 2), (1, 2)), 0).hsp.context,
                context
            );
        }
        let mut method = hsp(1, 1, (1, 2), (1, 2));
        method.comp_adjustment_method = 1;
        assert_eq!(blastx_hsp(&method, 0).composition_method, 1);
    }

    // A query's hit list: one BLASTX list per profile with HSPs, in the hit list's order, named
    // by the profile's ordinal id; the best e-value refreshed (`HspList::refresh`).
    #[test]
    fn a_hit_list_becomes_blastx_lists_in_its_order() {
        let list = |oid: u32, evalues: &[f64]| BlastpHspList {
            oid,
            query_index: 0,
            hsps: evalues
                .iter()
                .map(|&e| BlastpHsp {
                    e_value: e,
                    ..hsp(1, 10, (1, 9), (1, 9))
                })
                .collect(),
            best_evalue: 0.0,
        };
        let hit_list = BlastpHitList {
            hsplist_count: 3,
            hsplist_max: 3,
            worst_evalue: 0.0,
            low_score: 0,
            heapified: false,
            hsplist_array: vec![
                list(5, &[1e-3, 1e-9]),
                list(1, &[]),
                list(2, &[0.5]),
                list(9, &[1e-30]),
            ],
            hsplist_current: 0,
            num_hits: 0,
        };
        let lists = blastx_lists(Some(hit_list), 3);
        let oids: Vec<usize> = lists.iter().map(|l| l.oid).collect();
        assert_eq!(oids, [5, 2]);
        assert!(lists.iter().all(|l| l.query_index == 3));
        assert_eq!(lists[0].best_evalue, 1e-9);
        assert_eq!(lists[0].hsps.len(), 2);
        assert_eq!(lists[0].hsps[0].hsp.context, 3 * 6);
        assert!(blastx_lists(None, 0).is_empty());
    }

    // local_blast.cpp:176-181 (quoted at `search_batch`): a batch without a valid frame is not
    // searched; its results are BLASTX's skipped batch.
    #[test]
    fn a_batch_without_a_valid_frame_is_not_searched() {
        let options = options(&[]);
        let records = [record("n", b"NNNNNN"), record("a", b"AC")];
        let prepared = prepare_unit_queries(&records, None, &options).unwrap();
        let parameters = score_block(&prepared);
        assert!(!parameters.iter().any(|p| p.valid));
        let results = unsearched_batch(prepared, parameters);
        assert!(results.search_skipped);
        assert_eq!(results.queries.len(), 2);
        assert!(results.queries.iter().all(Vec::is_empty));
        assert_eq!(results.parameters.len(), 12);
    }

    // The refusals where the search starts.
    #[test]
    fn sum_statistics_and_an_ungapped_search_are_refused_where_the_search_starts() {
        let mut search = TranslatedRpsSearch {
            search: RpsSearch::default(),
            formatting: FormattingArgs::default(),
            ungapped: false,
        };
        assert!(search.check_search().is_ok());
        search.search.sum_stats = true;
        assert!(search
            .check_search()
            .unwrap_err()
            .to_string()
            .contains("-sum_stats true"));
        search.search.sum_stats = false;
        search.ungapped = true;
        assert!(search
            .check_search()
            .unwrap_err()
            .to_string()
            .contains("crashes after an ungapped search"));
        assert_eq!(
            search.formatting(),
            FormattingArgs::default(),
            "the formatting options reach BLASTX's reports"
        );
    }
}
