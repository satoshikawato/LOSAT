//! BLASTP's split query batches and NCBI's steps of a split batch (`SplitBatchSteps`, the one owner
//! of BLASTP's and TBLASTN's): which searches split their batches, a chunk's parts and validity,
//! and where the failures of the setup of a call's batches come when a split batch meets them first
//! (`first_batch_setup` for `run_call`, psiblast's first read, the set-up of a refused batch; line
//! RG item X-bpgap).

use anyhow::Result;

use super::call_setup::{self, CallSetup};
use super::database::{BlastpDatabase, PreliminaryReadError};
use super::{blastp_query_frames, pssm_hooks, BlastpSubjectSource, ResolvedBlastpArgs};
use crate::algorithm::blastn::query_split::{ChunkQuery, SplitSizes};
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::protein_query_split::{
    protein_chunk_part_in_record, protein_query_split_sizes, split_batch_reaches_read,
    split_protein_batch, splits_query_batches, ProteinQueryChunk, ReadFails, SetupCheckFailure,
    SplitBatchSteps, UnmetChecks,
};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::blastinput::app::ErrorAfterBatchReports;
use crate::blastinput::app::TRACEBACK_FAILED;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::query_batch::BatchSetUp;
use crate::blastinput::seq_range;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:86-93
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:61
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:122-129 (cited at
// `psi_round::run_queries`: a query factory that is not `CObjMgr_QueryFactory` makes one chunk)
/// Whether the searches of `args` split their query batches (`search_split_query_batches`), by the
/// one rule of BLASTP's, psiblast's and TBLASTN's (`splits_query_batches`): an RPS-BLAST query
/// batch is never split; an ungapped search splits none (and so builds one compressed lookup table
/// for each batch, `ncbi_compressed_lookup_tables`); a PSSM search splits only the rounds of a
/// FASTA query (`PsiSearch::split_setup`), never a PSSM query (`-in_pssm`,
/// `CObjMgrFree_QueryFactory`); PHI-BLAST (`-phi_pattern`) splits none.
pub(super) fn splits_batches(args: &ResolvedBlastpArgs) -> bool {
    splits_query_batches(
        !args.ungapped,
        args.psi_search
            .as_deref()
            .is_some_and(|psi| psi.split_setup.is_none()),
        args.rps_search.is_some(),
        args.phi.is_some(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:278-281
/// The record of the part `part` of the query `record` as its chunk reads it: with
/// `-lcase_masking`, the lower case of the query's whole record restricted to the part as NCBI
/// restricts it (`protein_chunk_part_in_record`: `whole`, the record and where the query's
/// searched letters start in it; each run one letter longer, and with `-query_loc` in the part's
/// coordinates without the interval's start).
pub(super) fn chunk_part(
    args: &ResolvedBlastpArgs,
    record: &FastaRecord,
    whole: (&[u8], usize),
    part: &ChunkQuery,
) -> FastaRecord {
    let mut cut = record.cut(part.from, part.to);
    if args.lcase_masking {
        cut.sequence = protein_chunk_part_in_record(whole.0, whole.1, part);
    }
    cut
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:191-196
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:273-284 (cited at
// `SplitBatchSteps`)
/// Whether a context of the query chunk of the parts `parts` is valid, as the chunk's setup finds
/// it: set up without the PSSM (a split round's chunks are searched as round 1's), with
/// `call_setup::call_lookup`'s validity (a context without letters is invalid). A chunk without
/// one finds nothing, as NCBI ignores its setup's error.
pub(super) fn chunk_valid(args: &ResolvedBlastpArgs, parts: &[FastaRecord]) -> bool {
    let mut unpssm;
    let args = match args.psi_search {
        Some(_) => {
            unpssm = args.clone();
            unpssm.psi_search = None;
            &unpssm
        }
        None => args,
    };
    pssm_hooks::psi_contexts(args, &blastp_query_frames(args, parts))
        .iter()
        .any(|context| context.is_valid && context.aa_len > 0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:206-215
/// Whether the search of a query chunk has a subject in hand (`SplitBatchSteps::before_search`):
/// a subject of the FASTA set with letters (one without is not searched, with its warning), or a
/// sequence that the database's lists include.
pub(super) fn subjects_in_hand(subject_source: BlastpSubjectSource<'_>) -> bool {
    match subject_source {
        BlastpSubjectSource::Fasta(records) => {
            records.iter().any(|record| !record.seq().is_empty())
        }
        BlastpSubjectSource::Database(database) => database.num_sequences() > 0,
    }
}

/// The form of a stop of the query batch whose first searched query is `query`: NCBI stops in the
/// search of the batch after it has printed the reports of the batches before (blastp_app.cpp:
/// 257-276, cited in `search`), the batch's failure (`ErrorAfterBatchReports::of_batch`, which
/// names the batch for `search`). psiblast searches one query at a time
/// and NCBI stops in the query's round 1 (it splits that round too), after the rounds of the
/// queries before: every round of a query makes the same chunks, so a later round never reaches a
/// stop, and the stop of the query's search is the plain error.
pub(super) fn batch_failure(
    args: &ResolvedBlastpArgs,
    query: usize,
    error: anyhow::Error,
) -> anyhow::Error {
    if args.psi_search.is_some() {
        error
    } else {
        // Before the batch's search: nothing of it read a header.
        ErrorAfterBatchReports::of_batch(query, error).into()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
/// The first query of the call's first searched batch (`batches`, with whether NCBI skips each
/// query's batch, `skipped`: no valid context), whose setup's failure (`batch_failure`) comes
/// after the reports of the batches before; 0 where the call searches no batch.
pub(super) fn first_searched_query(batches: &[seq_range::RangedBatch], skipped: &[bool]) -> usize {
    batches
        .iter()
        .find(|batch| !batch.searched.is_empty() && !skipped[batch.searched.start])
        .map_or(0, |batch| batch.searched.start)
}

/// The text of NCBI's engine error of a failed setup check of a batch (`SplitBatchSteps`' setup
/// checks), apart from LOSAT's refusals there.
#[derive(Debug)]
struct SetupText(&'static str);

impl std::fmt::Display for SetupText {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.0)
    }
}

impl std::error::Error for SetupText {}

/// The error that a failure of the setup of a batch `error` makes (a check's `SetupText`, or a read
/// that fails, `database::PreliminaryReadError`): NCBI's engine error; another error (LOSAT's
/// refusal of a count of 0, a reader's own error) as it is.
fn setup_failure(error: anyhow::Error) -> anyhow::Error {
    let error = match error.downcast::<SetupText>() {
        Ok(text) => return crate::blastinput::app::engine_error(text.0),
        Err(error) => error,
    };
    match error.downcast::<PreliminaryReadError>() {
        Ok(read) => crate::blastinput::app::engine_error(read.as_bytes()),
        Err(error) => error,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-982 (cited at
// `SplitBatchSteps`)
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:228-250
/// What the subjects give the setup of a batch: the read of OID 0 (`first_oid_read_fails`) or
/// the search's read of a sequence file that cannot be opened (`sequence_file_read_error`, which
/// `CLocalBlast` reports as `Database memory map file error`), and the checks of the setup's
/// subject length after LOSAT's refusal of a count of 0 (`call_setup::setup_fails`, with the
/// Gumbel block of `has_gumbel`).
struct SetupFailures {
    read_fails: Option<ReadFails>,
    read_error: Option<anyhow::Error>,
    checks: Result<()>,
}

impl SetupFailures {
    fn of(
        args: &ResolvedBlastpArgs,
        database: Option<&BlastpDatabase>,
        subjects: &[EncodedProtein],
        any_valid: bool,
    ) -> Self {
        if database.is_some_and(BlastpDatabase::first_oid_read_fails) {
            return Self {
                read_fails: Some(ReadFails::Setup),
                read_error: Some(
                    PreliminaryReadError::from_static("Database memory map file error").into(),
                ),
                checks: Ok(()),
            };
        }
        let read_error =
            database.and_then(|database| database.sequence_file_read_error(args.app.name()));
        let checks = call_setup::db_extent(args, database, subjects).and_then(|extent| {
            match call_setup::setup_fails(args, &extent, any_valid)? {
                true => Err(SetupText(crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT).into()),
                false => Ok(()),
            }
        });
        Self {
            read_fails: read_error.is_some().then_some(ReadFails::Search),
            read_error,
            checks,
        }
    }

    fn pending(&self) -> bool {
        self.read_fails.is_some() || self.checks.is_err()
    }

    /// The batch's setup checks for its `SplitBatchSteps` (the error kept here for a batch after):
    /// NCBI's engine error of the setup (`SetupText`) or LOSAT's refusal.
    fn checks(&self) -> Result<(), SetupCheckFailure> {
        match &self.checks {
            Ok(()) => Ok(()),
            Err(error) => Err(match error.downcast_ref::<SetupText>() {
                Some(text) => SetupCheckFailure::Engine(SetupText(text.0).into()),
                None => SetupCheckFailure::Refusal(anyhow::anyhow!("{error:#}")),
            }),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-982 (cited at
    // `SplitBatchSteps`)
    /// The failure that a batch searched whole with a valid context meets, in NCBI's order: the
    /// read of OID 0, the checks, the search's read.
    fn of_whole_batch(self) -> Option<anyhow::Error> {
        match (self.read_fails, self.read_error, self.checks) {
            (Some(ReadFails::Setup), Some(error), _) => Some(error),
            (_, _, Err(error)) => Some(error),
            (_, error, Ok(())) => error,
        }
    }
}

/// What the setup of the call's first searched batch does before its checks
/// (`first_batch_setup`) where a split batch is the call's first searched one (`call_steps`).
enum CallSteps {
    /// The call's own checks and search: its first searched batch meets the setup's failures as
    /// `first_batch_setup` finds them, if any; `reads`: whether a batch's search reads the database (none
    /// does where every searched batch is split into chunks without a valid context: the
    /// preliminary stage then reads nothing, as for a call of invalid queries).
    Search { reads: bool },
    /// NCBI's engine error of the call's first searched batch, split, whose chunks found nothing
    /// (`TRACEBACK_FAILED`, `preliminary_search_error`).
    FirstBatchError(&'static str),
    /// A later batch's failure or stop after the reports of the batches before
    /// (`batch_failure`), or LOSAT's refusal: the call's error.
    Fails(anyhow::Error),
}

impl CallSteps {
    /// Whether the call's search reads the database (`Search`), or the call's error, the first
    /// searched batch's engine error made of its text by `first_batch_error`.
    fn reads(self, first_batch_error: impl FnOnce(&'static str) -> anyhow::Error) -> Result<bool> {
        match self {
            CallSteps::Search { reads } => Ok(reads),
            CallSteps::FirstBatchError(text) => Err(first_batch_error(text)),
            CallSteps::Fails(error) => Err(error),
        }
    }
}

/// The batches of a call (`call_setup::call_lookup`): each batch, whether NCBI skips it (no valid
/// context, local_blast.cpp:177-180), and the queries' lengths.
#[derive(Clone, Copy)]
pub(super) struct CallBatches<'c>(
    pub(super) &'c [seq_range::RangedBatch],
    pub(super) &'c [bool],
    pub(super) &'c [usize],
);

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-982 (cited at
// `SplitBatchSteps`)
/// Where NCBI's searches of the call's batches meet the failures of their setup (the read of OID
/// 0, the checks of the subject length, the search's first read) when the call's first searched
/// batch is split (NCBI 2.17.0 probes of line RG item X-bpgap): the batches in order, a batch of
/// invalid queries only not searched, one searched whole meeting them as `first_batch_setup` does,
/// a split one meeting its chunks' stops first and the failures where a chunk gets to them
/// (`SplitBatchSteps`), or, where no chunk has a valid context, the traceback's own setup meeting
/// the checks (`Traceback failed`) or nothing read, the next batch going on. A call whose first
/// searched batch is searched whole, or whose setup does not fail, meets `first_batch_setup`'s
/// checks as before: a split batch's stops then come in `search_split_query_batches`, after the
/// call's search.
fn call_steps(
    call: &CallSetup<'_, '_>,
    subject_source: BlastpSubjectSource<'_>,
    subjects: &[EncodedProtein],
    contexts: &[QueryContext],
    batches: CallBatches<'_>,
) -> Result<CallSteps> {
    let args = call.args;
    let searches = CallSteps::Search { reads: true };
    if call.prelim.is_some() || !splits_batches(args) {
        return Ok(searches);
    }
    let CallBatches(query_batches, skipped, query_lengths) = batches;
    let searched: Vec<std::ops::Range<usize>> = query_batches
        .iter()
        .map(|batch| batch.searched.clone())
        .filter(|searched| !searched.is_empty() && !skipped[searched.start])
        .collect();
    let sizes = protein_query_split_sizes(10000)?;
    let split_of =
        |range: &std::ops::Range<usize>| split_protein_batch(&query_lengths[range.clone()], sizes);
    if searched.first().and_then(&split_of).is_none() {
        return Ok(searches);
    }
    let failures = SetupFailures::of(
        args,
        subject_source.database(),
        subjects,
        contexts.iter().any(|ctx| ctx.is_valid),
    );
    if !failures.pending() {
        return Ok(searches);
    }
    // The batches after the first come after split batches whose chunks had no valid context and
    // met no check (`searched_before`): the checks, the call's, passed, and such a batch meets the
    // reads that fail after the reports of the batches before.
    let mut searched_before = false;
    for range in searched {
        let after_before = |error: anyhow::Error| -> Result<CallSteps> {
            Ok(CallSteps::Fails(batch_failure(
                args,
                range.start,
                setup_failure(error),
            )))
        };
        let Some(chunks) = split_of(&range) else {
            // Not the first searched batch, which is split (above).
            return match failures.of_whole_batch() {
                Some(error) => after_before(error),
                None => Ok(searches),
            };
        };
        let stop = |refusal: anyhow::Error| batch_failure(args, range.start, refusal);
        let mut steps = match SplitBatchSteps::new(
            args.app.name(),
            &chunks,
            sizes,
            failures.read_fails,
            failures.checks(),
            &stop,
        ) {
            Ok(steps) => steps,
            Err(error) => return Ok(CallSteps::Fails(error)),
        };
        let parts = |chunk: &ProteinQueryChunk| {
            chunk_parts(
                args,
                (call.query_records, call.ranges.input_records),
                call.ranges.query_input,
                range.start,
                chunk,
            )
        };
        let reaches = split_batch_reaches_read(&mut steps, &chunks, |chunk| {
            Ok(chunk_valid(args, &parts(chunk)))
        });
        match reaches {
            // A stop of the chunks (in its form).
            Err(error)
                if error.is::<ErrorAfterBatchReports>()
                    || (args.psi_search.is_some() && !error.is::<SetupText>()) =>
            {
                return Ok(CallSteps::Fails(error))
            }
            // The checks of the first valid chunk of the first searched batch, which
            // `first_batch_setup` meets.
            Err(_) => {
                debug_assert!(!searched_before, "a later batch meets no check");
                return Ok(searches);
            }
            Ok(true) if searched_before => {
                return match failures.of_whole_batch() {
                    Some(error) => after_before(error),
                    None => Ok(searches),
                }
            }
            Ok(true) => return Ok(searches),
            Ok(false) => {}
        }
        match steps.unmet_checks() {
            Some(unmet) => {
                debug_assert!(!searched_before, "a later batch meets no check");
                // The refusal of a count of 0 that the traceback's own setup of this batch meets
                // (NCBI divides by the count there and crashes) comes after the reports of the
                // batches before, of invalid queries only (`batch_failure`, line RG item bp-SF1).
                return Ok(match unmet {
                    UnmetChecks::TracebackFailed => CallSteps::FirstBatchError(TRACEBACK_FAILED),
                    UnmetChecks::Refusal(refusal) => {
                        CallSteps::Fails(batch_failure(args, range.start, refusal))
                    }
                });
            }
            None => searched_before = true,
        }
    }
    Ok(CallSteps::Search {
        reads: failures.read_fails.is_none(),
    })
}

/// The setup of the call's first searched batch before the call's search (`run_call`), in NCBI's
/// order: a split first searched batch's steps (`call_steps`), then the read of OID 0 of a database
/// without sequences, the database's extent (`call_setup::db_extent`) and NCBI's error of its
/// setup's average subject length (`call_setup::setup_fails`, after LOSAT's refusal of a count of
/// 0, the first searched batch's failure: `batch_failure`). `setup_error` makes NCBI's engine
/// error of the first searched batch
/// (`preliminary_search_error`) of a database's name and a text. Gives whether the call's search
/// reads the database (`CallSteps::Search`) and the extent.
pub(super) fn first_batch_setup<'a>(
    call: &CallSetup<'_, '_>,
    (subject_source, subjects): (BlastpSubjectSource<'a>, &[EncodedProtein]),
    contexts: &[QueryContext],
    batches: CallBatches<'_>,
    mut setup_error: impl FnMut(&[u8], &str) -> anyhow::Error,
) -> Result<(bool, call_setup::BlastpDbExtent<'a>)> {
    let (args, database, skipped) = (call.args, subject_source.database(), batches.1);
    let reads = call_steps(call, subject_source, subjects, contexts, batches)?
        .reads(|text| setup_error(b"", text))?;
    // A database without sequences (makeblastdb builds none; a volume whose index holds
    // no sequence) and without an alias file's LENGTH or STATS_TOTLEN
    // (`SeqDb::first_oid_read_fails`): NCBI opens it, writes the outfmt 0 prolog, reads the
    // first query batch and fails in that batch's preliminary search, where its reader cannot
    // map the sequence file (fixture `LOSAT/tests/seqdb/db/empty`).
    // A first batch of invalid queries only is not searched: the batches up to the first one that
    // is are reported before the error (`unsearched_batches_then_engine_error`); a call of invalid
    // queries only searches no batch, and NCBI reports them (NCBI 2.17.0 blastp, an all-X query
    // and the empty database: exit 0).
    if let Some(database) = database.filter(|database| {
        database.first_oid_read_fails()
            && (call.query_chunk || skipped.iter().any(|skipped| !skipped))
    }) {
        return Err(setup_error(
            database.name(),
            "Database memory map file error",
        ));
    }
    let extent = call_setup::db_extent(args, database, subjects)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276
    // LOSAT's refusal of a count of 0 (`call_setup::setup_fails`) is the failure of NCBI's first
    // searched batch, whose setup divides by the count and crashes (SIGFPE) after the batch loop
    // has printed the reports of the batches before, of invalid queries only: the refusal comes
    // after them (`batch_failure`; line RG item bp-SF1, NCBI 2.17.0 probes of six X-only batches
    // before a valid query, as BLASTN's and TBLASTN's of D1a-5).
    let any_valid = contexts.iter().any(|ctx| ctx.is_valid);
    let fails = call_setup::setup_fails(args, &extent, any_valid).map_err(|refusal| {
        batch_failure(args, first_searched_query(batches.0, skipped), refusal)
    })?;
    // The setup of the first searched batch fails: NCBI has written the outfmt 0 prolog and the
    // batch's reader messages (`preliminary_search_error`; NCBI 2.17.0 blastp and psiblast
    // `-subject` of the letters `M`, ``, `` with `-searchsp` or `OLD_FSC`: exit 3). After a first
    // batch of invalid queries only, which NCBI does not search, the error follows the reports of
    // the batches before the first searched one (`unsearched_batches_then_engine_error`).
    if fails {
        return Err(setup_error(
            b"",
            crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT,
        ));
    }
    Ok((reads, extent))
}

/// The parts of the query chunk `chunk` of the batch whose first searched query is `batch_start`,
/// of the queries `query_records` (`chunk_part`; their whole records `input_records`, by
/// `input`'s ordinals and placements).
pub(super) fn chunk_parts(
    args: &ResolvedBlastpArgs,
    (query_records, input_records): (&[FastaRecord], &[FastaRecord]),
    input: &seq_range::QueryInput,
    batch_start: usize,
    chunk: &ProteinQueryChunk,
) -> Vec<FastaRecord> {
    chunk
        .queries
        .iter()
        .map(|part| {
            let query = batch_start + part.query;
            let whole = (
                input_records[input.ordinal(query)].seq(),
                input.placements.offset(query),
            );
            chunk_part(args, &query_records[query], whole, part)
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:233-290 (cited at
// `SplitBatchSteps`)
/// psiblast: whether the search of the query `record` (with its whole record `whole`), which NCBI
/// splits, gets to the first read of the database, whose read fails in `read_fails`
/// (`psiblast_stop_at_invalid_query`), or its chunks' stop (the plain error); `None` for a query
/// that NCBI does not split.
pub(super) fn split_query_reads(
    args: &ResolvedBlastpArgs,
    record: &FastaRecord,
    whole: (&[u8], usize),
    read_fails: ReadFails,
) -> Option<Result<bool>> {
    if !splits_batches(args) {
        return None;
    }
    let sizes = match protein_query_split_sizes(10000) {
        Ok(sizes) => sizes,
        Err(error) => return Some(Err(error)),
    };
    let chunks = split_protein_batch(&[record.seq().len()], sizes)?;
    Some(reaches_read(
        args,
        &chunks,
        sizes,
        Some(read_fails),
        |chunk| {
            let parts: Vec<FastaRecord> = chunk
                .queries
                .iter()
                .map(|part| chunk_part(args, record, whole, part))
                .collect();
            chunk_valid(args, &parts)
        },
    ))
}

/// `split_batch_reaches_read` of the chunks `chunks` with plain stops and checks that pass.
fn reaches_read(
    args: &ResolvedBlastpArgs,
    chunks: &[ProteinQueryChunk],
    sizes: SplitSizes,
    read_fails: Option<ReadFails>,
    mut valid: impl FnMut(&ProteinQueryChunk) -> bool,
) -> Result<bool> {
    let plain = |refusal: anyhow::Error| refusal;
    let mut steps =
        SplitBatchSteps::new(args.app.name(), chunks, sizes, read_fails, Ok(()), &plain)?;
    split_batch_reaches_read(&mut steps, chunks, |chunk| Ok(valid(chunk)))
}

/// A chunk stop of the set-up of a refused batch (`refused_batch_setup_error`), which keeps the
/// refusal.
#[derive(Debug)]
struct ChunkStop;

impl std::fmt::Display for ChunkStop {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("a stop of a query chunk")
    }
}

impl std::error::Error for ChunkStop {}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-982 (cited at
// `SplitBatchSteps`)
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:233-290 (cited at
// `SplitBatchSteps`)
/// BLASTP's set-up of a refused batch (`blast_engine::refused_batch_stop`; TBLASTN's
/// `stage_d_pipeline::split_batch_reads`, review 4 of line RG item X-short-bat, M1): NCBI's batch
/// is `records` (each searched query's record and, where it is first in its scope, its whole
/// record, the records set up). NCBI's error where a context of theirs is valid and the set-up
/// fails for the subjects (`SetupFailures`, as `run_call` fails the first searched batch: a
/// database without sequences, `Database memory map file error`; the average subject length
/// without the Gumbel block, `-searchsp`, `OLD_FSC` or an ungapped search), or where the search's
/// first read of the database then fails (a sequence file that cannot be opened, `Database memory
/// map file error`; line WQA D10 round 6, review 4 M2), or LOSAT's (a count of 0), which keeps the
/// batch's refusal (`query_batch::refused_batch_error`); `Searched` where it passes with a valid
/// context, so that the refusal comes unless the run's XML failure for its matrix name takes it;
/// `Unknown` where no context is valid. A batch that NCBI splits (`splits_batches`) meets them in
/// its chunks' steps (NCBI 2.17.0 probes `d1`, `d2` of line RG item X-bpgap), where a chunk's stop
/// keeps the refusal (`Unknown`): a chunk is valid where a part of a record first in its scope has
/// a valid context (the parts of a repeated id, which NCBI reads in its earlier record's letters,
/// are left out), and a batch whose chunks have none ends at the traceback's own setup (`Traceback
/// failed`).
pub(super) fn refused_batch_setup_error(
    args: &ResolvedBlastpArgs,
    database: Option<&BlastpDatabase>,
    subject_records: &[FastaRecord],
    records: &[(&FastaRecord, Option<(&[u8], usize)>)],
) -> BatchSetUp {
    let Ok(sizes) = protein_query_split_sizes(10000) else {
        return BatchSetUp::Unknown;
    };
    let lengths: Vec<usize> = records
        .iter()
        .map(|(record, _)| record.seq().len())
        .collect();
    let chunks = splits_batches(args)
        .then(|| split_protein_batch(&lengths, sizes))
        .flatten();
    let encoded: Vec<EncodedProtein> = match database {
        Some(_) => Vec::new(),
        None => subject_records
            .iter()
            .map(|record| crate::algorithm::blastp::encoding::encode_protein_sequence(record.seq()))
            .collect(),
    };
    let own: Vec<FastaRecord> = records
        .iter()
        .filter(|(_, whole)| whole.is_some())
        .map(|&(record, _)| record.clone())
        .collect();
    // `call_setup::call_lookup`'s validity: a context without letters is invalid.
    let contexts = pssm_hooks::psi_contexts(args, &blastp_query_frames(args, &own));
    if !contexts
        .iter()
        .any(|context| context.is_valid && context.aa_len > 0)
    {
        if chunks.is_some() {
            return BatchSetUp::Unknown;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-982
        // The subjects' failures of a batch with a valid context (`SetupFailures`, which read the
        // subjects alone): where there is none, NCBI's batch fails the same whatever the letters
        // that NCBI reads for its repeats (`BatchSetUp::NoValidContext`).
        let subjects_pass = SetupFailures::of(args, database, &encoded, true)
            .of_whole_batch()
            .is_none();
        return BatchSetUp::NoValidContext { subjects_pass };
    }
    let failures = SetupFailures::of(args, database, &encoded, true);
    let Some(chunks) = chunks else {
        return failures
            .of_whole_batch()
            .map_or(BatchSetUp::Searched, |failure| {
                BatchSetUp::Fails(setup_failure(failure))
            });
    };
    let stop = |refusal: anyhow::Error| refusal.context(ChunkStop);
    let Ok(mut steps) = SplitBatchSteps::new(
        args.app.name(),
        &chunks,
        sizes,
        failures.read_fails,
        failures.checks(),
        &stop,
    ) else {
        return BatchSetUp::Unknown;
    };
    let reaches = split_batch_reaches_read(&mut steps, &chunks, |chunk| {
        let parts: Vec<FastaRecord> = chunk
            .queries
            .iter()
            .filter_map(|part| {
                let (record, whole) = records[part.query];
                whole.map(|whole| chunk_part(args, record, whole, part))
            })
            .collect();
        Ok(!parts.is_empty() && chunk_valid(args, &parts))
    });
    match reaches {
        Err(error) if error.is::<ChunkStop>() => BatchSetUp::Unknown,
        Err(error) => BatchSetUp::Fails(setup_failure(error)),
        Ok(true) => failures
            .of_whole_batch()
            .map_or(BatchSetUp::Searched, |failure| {
                BatchSetUp::Fails(setup_failure(failure))
            }),
        Ok(false) => match steps.unmet_checks() {
            Some(UnmetChecks::TracebackFailed) => {
                BatchSetUp::Fails(crate::blastinput::app::engine_error(TRACEBACK_FAILED))
            }
            Some(UnmetChecks::Refusal(refusal)) => BatchSetUp::Fails(refusal),
            None => BatchSetUp::Unknown,
        },
    }
}
