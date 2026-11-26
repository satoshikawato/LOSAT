//! The stages of an rpsblast search (series PR, item PS11a) that differ from blastp's, called
//! once per run from `search_cli` and `run_call` when the search is RPS-BLAST's
//! (`ResolvedBlastpArgs::rps_search`): the check of the database's frequency ratios, the queries'
//! contexts, and the whole RPS search (`algorithm::rps_search`) into the final hit lists, whose
//! subjects are the profiles' consensus sequences of the database. blastp's driver does the rest
//! (the query batches, the database's statistics, the reports). Every function is out of line
//! and cold: blastp's default path only tests `rps_search`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:129-141

use crate::wording::ENGINE_ERROR;
use anyhow::{bail, Result};

use super::{database, BlastpRuntimeLookup};
use crate::algorithm::blastp::args::{BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::rps_search::engine::{
    batch_gapped_k, RpsFrame, RpsQuery, RpsRun, RpsSearch, RpsWorkers,
};
use crate::algorithm::rps_search::volumes::{
    mask_unsupported_letters, merge_volume, report_ordinal_ids, shift_ordinal_ids, ProfileLetters,
    RpsUnits, RpsVolumeError, RpsVolumeSet, RpsVolumeSetup,
};
use crate::algorithm::tblastx::lookup::{
    prepare_blosum62_lookup_query_for_word_size, QueryContext,
};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::stats::KarlinParams;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:158-164
/// The error of the RPS files that NCBI reports itself (`BLAST engine error: <text>`, exit
/// status 3): a volume's files that cannot be mapped or read, or the score block of an `.aux`
/// without a matrix name (review B of PS12a, N4). NCBI sets up every volume's files in the first
/// batch's set-up (`BlastSetupPreliminarySearchEx`, before `CheckInternalData` looks at the
/// batch's queries), so the error is the first batch's whatever its queries' validity, after
/// what that batch writes before (`preliminary_search_error`; series PR, item PS12f, review S1).
#[derive(Debug)]
pub(super) struct RpsEngineError(pub(super) Vec<u8>);

impl std::fmt::Display for RpsEngineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{ENGINE_ERROR} {}", String::from_utf8_lossy(&self.0))
    }
}

impl std::error::Error for RpsEngineError {}

/// `CheckForFreqRatioFile`: with composition-based statistics, each database name of `-db` whose
/// first volume has no `.freq` file gets NCBI's warning and the search runs without them
/// (`args.comp_based_stats` becomes mode 0). Only the first volume of each name is looked at: a
/// later volume without `.freq` fails inside the search (`RpsVolumeFiles::open_freq`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:807-840
#[cold]
#[inline(never)]
pub(super) fn check_freq_ratio_file(
    args: &mut ResolvedBlastpArgs,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    diagnostics: &mut dyn std::io::Write,
) -> Result<()> {
    if !args.comp_based_stats.is_enabled() {
        return Ok(());
    }
    let Some(name) = args.db.as_deref() else {
        return Ok(());
    };
    let path = crate::seqdb::SearchPath::from_environment(settings.blastdb.as_deref());
    let mut failed: Vec<&str> = Vec::new();
    for db in name.split(' ').filter(|db| !db.is_empty()) {
        let volumes = crate::seqdb::SeqDb::volume_and_alias_paths(
            db.as_bytes(),
            crate::seqdb::SeqType::Protein,
            &path,
        )
        .ok()
        .map(|(volumes, _)| volumes)
        .filter(|volumes| !volumes.is_empty())
        .ok_or_else(|| {
            crate::blastinput::app::engine_error("Cannot retrieve path to RPS database")
        })?;
        let freq = [volumes[0].as_slice(), b".freq"].concat();
        if !crate::utils::os_bytes::path_of_bytes(&freq).exists() {
            failed.push(db);
        }
    }
    if failed.is_empty() {
        return Ok(());
    }
    args.comp_based_stats.mode =
        crate::algorithm::blastp::args::BlastpCompositionMode::NoCompositionBasedStats;
    let warning = format!(
        "Warning: [rpsblast] {} contain(s) no freq ratios needed for composition-based statistics.\nRPSBLAST will be run without composition-based statistics.\n",
        failed.join(", ")
    );
    match &args.archive {
        Some(handle) => handle.write_warning(diagnostics, warning.as_bytes())?,
        None => diagnostics.write_all(warning.as_bytes())?,
    }
    Ok(())
}

/// What rpsblast sets in the options handle that the requests of the archive and of the search
/// strategy hold
/// (`report::archive::OptionsHandle::set_search_sets`): before the batches, the composition-based
/// statistics that `check_freq_ratio_file` turned off (NCBI's `SetCompositionBasedStats` in place:
/// `CompositionBasedStats 0` in every request, NCBI BLAST+ 2.17.0 `p9.rps.err.nofreq`); at each
/// batch's set-up, the database's length and `.aux` costs (`rps_search::volumes::search_sets`).
#[cold]
#[inline(never)]
pub(super) fn note_handle_sets(
    handle: &crate::report::archive::OptionsHandle,
    rps: &RpsSearch,
    args: &ResolvedBlastpArgs,
    database: &database::BlastpDatabase,
) -> Result<()> {
    use crate::report::archive::options::{Param, ParamValue};
    let before_batches = (rps.requested_cbs != 0 && !args.comp_based_stats.is_enabled())
        .then_some(Param {
            name: "CompositionBasedStats",
            value: ParamValue::Integer(0),
        })
        .into_iter()
        .collect();
    // The volumes' matrix and gap costs as the first volume's files read them (the search sets the
    // volumes up again, `rps_search`), after `x_AdjustDbSize`'s length (`RpsUnits`).
    let db_length = database.rps_units().and_then(RpsUnits::strategy_db_length);
    let volumes = RpsVolumeSet::of(database.seqdb());
    let costs = RpsVolumeSetup::new(rps.program()).costs(volumes.first());
    handle.set_search_sets(
        before_batches,
        crate::algorithm::rps_search::volumes::search_sets(db_length, costs),
    )
}

/// `database` with the search units of rpsblast's search of it and their statistics
/// (`rps_search::volumes::RpsUnits`, one owner with rpstblastn's, built once per run here, where
/// the database opens: `open_search_database`), with the run's `-searchsp` and `-dbsize`. Their
/// refusals come where NCBI meets them: the setup's counts in `run_call` where a query of the run
/// is valid (`call_setup::setup_fails`), the concatenated lengths at each unit's search
/// (`rps_search`, `RpsUnit::check_engine_length`).
#[cold]
#[inline(never)]
pub(super) fn with_units(
    args: &ResolvedBlastpArgs,
    database: database::BlastpDatabase,
) -> Result<database::BlastpDatabase> {
    let program = args
        .rps_search
        .as_deref()
        .map_or("RPSBLAST", RpsSearch::program);
    let units = RpsUnits::of(
        database.seqdb(),
        args.searchsp.unwrap_or(0),
        args.dbsize,
        program,
    )?;
    Ok(database.with_rps_units(units))
}

/// The queries' contexts (their letters and `kbp_std`, the composition with BLOSUM62, and
/// validity) and an empty lookup table: the RPS search scans with the database's prebuilt table.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2746-2748
#[cold]
#[inline(never)]
pub(super) fn rps_contexts(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    let (_, _, contexts) =
        prepare_blosum62_lookup_query_for_word_size(query_frames, args.word_size, false);
    (super::empty_aa_lookup(args.word_size), contexts)
}

/// The per-run values of the search that `run_call` computed for blastp.
pub(super) struct RpsSearchRun<'a> {
    pub(super) contexts: &'a [QueryContext],
    /// `kbp_gap_std`: BLOSUM62 11/1's table values.
    pub(super) gapped_params: &'a KarlinParams,
    /// Whether the search has the Gumbel block (none with `-searchsp`, at any value, or
    /// `OLD_FSC`): each unit's is made from its length (`RpsUnits::gumbel_length`).
    pub(super) gumbel: bool,
    /// The search units of the run and their statistics (`units`; `None` without a database):
    /// each query's search space, the Gumbel block and the hit cutoff's subject length in each
    /// volume's search.
    pub(super) units: Option<&'a RpsUnits>,
    /// Each query's effective search space in the first unit, which the set-up computed with that
    /// unit's numbers (`RpsSetup`; `exact_search_spaces`, `None` where the length is 0).
    pub(super) search_spaces: &'a [Option<i64>],
    pub(super) prelim_hitlist_size: usize,
    /// The query batches (searched-query ranges) of the run.
    pub(super) batches: Vec<std::ops::Range<usize>>,
    /// Whether the engine counts each HSP's identities again for the reports: NCBI's reports
    /// count them on the sequences of their scope (outfmt 0 always, the tabular formats inside
    /// their alignment block), not on the masked query of the traceback's count
    /// (`rps_search::traceback::formatter_identities`). LOSAT's outfmt 0 and default tabular
    /// writers print the HSP's count, so the engine counts again for them; the tabular writer of
    /// a field list or outfmt 10 counts on the rows itself and reads the traceback's count for
    /// the block's gate (`tabular::needs_rows`). Decided once from the reports that the call
    /// writes (`call_setup::requested_reports`: no `BlastpReport::CustomTabular`), so that the
    /// count and the report read one decision.
    pub(super) recount_identities: bool,
}

/// What blastp's set-up of an rpsblast search reads of its first search unit in place of the
/// database's numbers (`run_call`): the merged result keeps that unit's statistics
/// (`thread_results[0]`), whose search spaces the footers print; each unit's search reads its own
/// (`rps_search`).
#[derive(Clone, Copy)]
pub(super) struct RpsSetup {
    /// `BLAST_GapAlignSetUp`'s total and count (`RpsUnit::setup_total`, `setup_count`): the
    /// average of the hit cutoff without the Gumbel block and the "too short" test.
    pub(super) setup_total: i64,
    pub(super) setup_count: i32,
    /// `BLAST_CalcEffLengths`' length and count (`RpsUnits::eff_lengths`).
    pub(super) eff_db_length: i64,
    pub(super) eff_num_seqs: i32,
    /// The Gumbel block's length (`RpsUnits::gumbel_length`, `-dbsize` included).
    pub(super) gumbel_db_length: i64,
    /// The hit cutoff's subject length with the Gumbel block (`RpsUnit::min_length`).
    pub(super) min_length: i32,
}

impl RpsSetup {
    /// The first unit of `units`.
    pub(super) fn of(units: &RpsUnits) -> Self {
        let first = units.first();
        let (eff_db_length, eff_num_seqs) = units.eff_lengths(first);
        Self {
            setup_total: first.setup_total,
            setup_count: first.setup_count,
            eff_db_length,
            eff_num_seqs,
            gumbel_db_length: units.gumbel_length(first),
            min_length: first.min_length,
        }
    }
}

/// The number of volumes that an RPS search of the database runs (`-num_threads 0` runs a thread
/// for each: NCBI's "auto"). A volume that the alias files reach twice is one volume of the
/// database's set (`SeqDb`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:414-424
#[cold]
#[inline(never)]
pub(super) fn volume_count(database: &database::BlastpDatabase) -> usize {
    database.seqdb().volumes().len().max(1)
}

/// What the RPS search of the run returns besides the hit lists.
pub(super) struct RpsSearched {
    /// The gapped K of each query's outfmt 0 footer.
    pub(super) gapped_k: Vec<f64>,
    /// The number of volumes searched.
    pub(super) volumes: usize,
    /// The HSPs of each query's profiles that the alias files' OID range leaves out past the
    /// formatter's cut (`report_ordinal_ids`), which left the hit lists: outfmt 7's count and the
    /// prefetch's limit count them (`report_prep::report_values`). Empty where none left.
    pub(super) unlisted: Vec<usize>,
}

/// The RPS search of every query of the run against each volume of the database in turn: the
/// final hit lists of the volumes merged into `hit_lists` (by query; their subjects named by the
/// database's ordinal ids, `report_ordinal_ids`; with several volumes longer than
/// `-max_target_seqs` and not to be sorted or cut by blastp's driver before the formatter's cut),
/// and the gapped K of each query's outfmt 0 footer (volume 0's).
///
/// Every volume is searched alone, the queries in parallel, and its files are unmapped before
/// the next volume's are mapped. NCBI searches the volumes in the sorted order of their paths
/// (a `set<string>`, `FindVolumePaths`), each a search unit with its own statistics
/// (`RpsUnits`: the search spaces, the Gumbel block and the shortest length under
/// `x_AdjustDbSize`'s length and count), a per-volume hit list size, and merges each query's
/// results into volume 0's (`s_MergeAlignSet`, `merge_volume`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:341-346
#[cold]
#[inline(never)]
pub(super) fn rps_search(
    rps: &RpsSearch,
    args: &ResolvedBlastpArgs,
    run: &RpsSearchRun<'_>,
    database: Option<&database::BlastpDatabase>,
    pool: &crate::utils::threading::SearchPool<'_>,
    hit_lists: &mut [Option<BlastpHitList>],
) -> Result<RpsSearched> {
    let program = rps.program();
    let (Some(database), Some(units)) = (database, run.units) else {
        bail!("a search without -db is not supported by LOSAT's {program}");
    };
    let db = database.seqdb();
    let RpsVolumeSet {
        volumes: db_volumes,
        starts,
        order,
    } = RpsVolumeSet::with_order(db, units.order().to_vec());
    let cbs = match args.comp_based_stats.mode {
        BlastpCompositionMode::NoCompositionBasedStats => 0u8,
        BlastpCompositionMode::CompositionBasedStats => 1,
        BlastpCompositionMode::CompositionMatrixAdjust => 2,
        BlastpCompositionMode::ForceFullMatrixAdjust => 3,
    };
    let searchsp = args.searchsp.unwrap_or(0);
    // Each query's effective search space in the unit being searched, with the effective lengths
    // of the unit it was computed for (`RpsUnits::eff_lengths`): the set-up's are the first unit's;
    // a unit of other lengths computes its own (`RpsUnits::search_space`; in the common case
    // every unit's lengths are the first's).
    let mut spaces: ((i64, i32), Vec<i64>) = (
        units.eff_lengths(units.first()),
        run.contexts
            .iter()
            .zip(run.search_spaces)
            .map(|(ctx, &space)| {
                if ctx.is_valid && ctx.aa_len > 0 {
                    space.unwrap_or(0)
                } else {
                    0
                }
            })
            .collect(),
    );
    let lambda = run.gapped_params.lambda;
    // BLOSUM62's score block for the matrix information of the redo.
    let sbp = if cbs > 0 {
        Some(crate::algorithm::rps_search::redo::redo_score_block(
            args.scoring.matrix,
            crate::core::psi::KarlinBlk {
                lambda,
                k: run.gapped_params.k,
                h: run.gapped_params.h,
            },
        )?)
    } else {
        None
    };
    // The X-drops in bits and the gap trigger (`RpsSearch::x_drop_bits`); the engine converts
    // them at each volume's scale (`RpsRun::cutoffs`, `RpsRun::gapped_x_drops`).
    let x_drops = rps.x_drop_bits();
    // The reports' identities (`RpsSearchRun::recount_identities`).
    let recount_identities = run.recount_identities;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250
    // A query that NCBI's engine would cut into chunks (`s_GetNextSubjectChunk`) is refused.
    if let Some(long) = run.contexts.iter().find(|ctx| ctx.aa_len > MAX_DBSEQ_LEN) {
        bail!(
            "a query of {} letters (NCBI BLAST+ searches a query above {MAX_DBSEQ_LEN} letters in chunks) is not supported by LOSAT's {program}",
            long.aa_len
        );
    }
    // The first query of each query's batch: its length (`-qcov_hsp_perc`) and, among the batch's
    // valid queries, the ungapped block that the redo's score block copy holds (`kbp_psi[0]`).
    let mut first_query_length = vec![0usize; run.contexts.len()];
    let mut batch_kbp: Vec<KarlinParams> =
        run.contexts.iter().map(|ctx| ctx.karlin_params).collect();
    for batch in &run.batches {
        if let Some(first) = run.contexts.get(batch.start) {
            for length in &mut first_query_length[batch.clone()] {
                *length = first.aa_len;
            }
        }
        if let Some(valid) = run.contexts[batch.clone()]
            .iter()
            .find(|ctx| ctx.is_valid && ctx.aa_len > 0)
        {
            for kbp in &mut batch_kbp[batch.clone()] {
                *kbp = valid.karlin_params;
            }
        }
    }
    let mut gapped_k = Vec::new();
    // The volumes' set-up kept from search to search of the database (`RpsRunState`); worker
    // scratches (`RpsWorkers`) live only for this search and are dropped before prefetch/reports
    // or DELTABLAST's main `-db` rounds.
    let mut state = database.rps_state();
    let state = &mut *state;
    let setup = state.setup.get_or_insert_with(|| {
        RpsVolumeSetup::new(program).with_database(database.name(), db_volumes.len())
    });
    let workers = RpsWorkers::new(pool.threads());
    // The first failing batch (its number and error) of the volumes searched so far.
    let spans: Vec<usize> = run.batches.iter().map(|batch| batch.len()).collect();
    let mut failed: Option<(usize, anyhow::Error)> = None;
    // Whether a query of the run is searched (NCBI's engine runs for a batch with a valid query).
    let searches = run
        .contexts
        .iter()
        .any(|ctx| ctx.is_valid && ctx.aa_len > 0);
    for ((position, &index), unit) in order.iter().enumerate().zip(units.units()) {
        let volume = &db_volumes[index];
        // The unit's Gumbel block (none with `-searchsp`, at any value, or `OLD_FSC`), from its
        // own length (`RpsUnits::gumbel_length`).
        let gumbel = if run.gumbel {
            Some(
                crate::stats::lookup_protein_gumbel_params(
                    &args.scoring,
                    units.gumbel_length(unit),
                )
                .ok_or_else(|| {
                    anyhow::anyhow!("missing NCBI gumbel parameters for supported blastp scoring")
                })?,
            )
        } else {
            None
        };
        let files = match setup.open_in_pool(volume, cbs > 0, Some(pool), rps.domain_search)? {
            Ok(files) => files,
            Err(RpsVolumeError::Files(text)) => return Err(RpsEngineError(text).into()),
            // N2 of review A (PS12a): the score block of each query of the first batch fails.
            Err(RpsVolumeError::NoMatrix) => {
                let queries = run.batches.first().map_or(0, |batch| batch.len());
                return Err(RpsEngineError(
                    "Error: Unknown error code -1 ".repeat(queries).into_bytes(),
                )
                .into());
            }
        };
        // The unit's concatenated length where NCBI's engine meets it, in the unit's search of the
        // first batch with a valid query (after the queries' set-up and the unit's setup, which
        // `run_call` checked): LOSAT's refusal (`RpsUnit::check_engine_length`; none where no
        // query is searched).
        if searches {
            unit.check_engine_length(args.window_size as i32, program)?;
        }
        let masked: Vec<Option<Vec<u8>>> = run
            .contexts
            .iter()
            .map(|ctx| mask_unsupported_letters(&ctx.aa_seq, files.alphabet))
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:932-940
        // (The stretch is NCBI's only inside the Gumbel block's branch: blast_parameters.c:922-945.)
        let prelim_evalue = if cbs > 1 && gumbel.is_some() {
            5.0 * args.evalue
        } else {
            args.evalue
        };
        // Each query's effective search space in this unit (`spaces`).
        let lengths = units.eff_lengths(unit);
        if lengths != spaces.0 {
            spaces = (
                lengths,
                run.contexts
                    .iter()
                    .map(|ctx| {
                        if !(ctx.is_valid && ctx.aa_len > 0) {
                            Ok(0)
                        } else {
                            units.search_space(
                                unit,
                                ctx.aa_len as i64,
                                run.gapped_params,
                                searchsp,
                                program,
                            )
                        }
                    })
                    .collect::<Result<Vec<i64>>>()?,
            );
        }
        let search_spaces = &spaces.1;
        // The preliminary stage's spaces where its lengths are not the setup's
        // (`RpsUnits::prelim_lengths`: one volume whose total is negative; `RpsRun::prelim_spaces`).
        let prelim_spaces = units
            .prelim_lengths(unit, searchsp)
            .map(|(length, count)| {
                run.contexts
                    .iter()
                    .map(|ctx| {
                        if !(ctx.is_valid && ctx.aa_len > 0) {
                            Ok(0)
                        } else {
                            Ok(crate::stats::setup::calc_eff_lengths(
                                ctx.aa_len as i64,
                                length,
                                count,
                                run.gapped_params,
                                0,
                                program,
                            )?
                            .eff_searchsp)
                        }
                    })
                    .collect::<Result<Vec<i64>>>()
            })
            .transpose()?;
        let engine = RpsRun {
            queries: run
                .contexts
                .iter()
                .zip(&first_query_length)
                .zip(&masked)
                .zip(&batch_kbp)
                .zip(search_spaces)
                .map(
                    |((((ctx, &first), masked), &batch_kbp), &search_space)| RpsQuery {
                        frames: vec![RpsFrame {
                            frame: 0,
                            search_space,
                            aa_seq: masked.as_deref().unwrap_or(&ctx.aa_seq),
                            valid: ctx.is_valid && ctx.aa_len > 0,
                            kbp_std: ctx.karlin_params,
                            first_query_length: first,
                            // The letters as read where the searched ones are masked: before
                            // SEG and the lower-case masks (`aa_seq_nomask`), else before the
                            // 26-column mask. The tabular report with a field list counts its
                            // identities on the rows itself, inside NCBI's alignment block,
                            // which the engine's count opens (`tabular::needs_rows`): it keeps
                            // that count.
                            aa_seq_nomask: if recount_identities {
                                ctx.aa_seq_nomask
                                    .as_deref()
                                    .or(masked.as_ref().map(|_| ctx.aa_seq.as_slice()))
                            } else {
                                None
                            },
                        }],
                        batch_kbp,
                    },
                )
                .collect(),
            gapped: *run.gapped_params,
            gumbel: gumbel.as_ref(),
            // The hit cutoff's subject length with the Gumbel block: the unit's shortest length
            // (without the block the cutoff reads the search space alone).
            cutoff_subject_length: unit.min_length,
            evalue: args.evalue,
            prelim_evalue,
            cbs,
            sbp: sbp.as_ref(),
            hitlist_size: args.traceback_hitlist_size(),
            prelim_hitlist_size: run.prelim_hitlist_size,
            window: args.window_size as i32,
            // The unit's concatenated length `Lc` (an alias file's `LENGTH` and `NSEQ` over one
            // volume).
            concat_length: unit.concat_length(),
            sorts_results: unit.sorts_results,
            prelim_spaces: prelim_spaces.as_deref(),
            x_drop_ungapped_bits: x_drops.ungapped,
            x_drop_gapped_bits: x_drops.gapped,
            x_drop_final_bits: x_drops.final_,
            gap_trigger_bits: x_drops.gap_trigger,
            options: rps,
            program,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:559-561
        // With sum statistics (`-sum_stats true`) NCBI's RPS search links the HSPs of each
        // profile with the user query's one-context `query_info` while the HSPs carry the
        // profile's context (0 to the number of profiles), so `BLAST_LinkHsps` reads past its
        // arrays (valgrind: invalid reads; a segmentation fault on some queries; e-values of 0.0
        // otherwise): undefined behaviour, refused where the first volume's search would start,
        // after NCBI's errors of the options, the database and its files.
        if rps.sum_stats {
            bail!(
                "-sum_stats true (NCBI BLAST+ 2.17.0's rpsblast links HSPs with the profiles' contexts as query indexes and reads past its arrays: its output depends on memory contents or it crashes) is not supported by LOSAT's {program}"
            );
        }
        // The run's query batches (`run.batches`) in NCBI's order: each batch's first refusal of
        // its cutoffs, else its first search error in query order (`RpsVolumeResults::batches`).
        // NCBI finishes batch 0 over every volume before it sets up batch 1, so a failing batch
        // ends the run only when it is batch 0; otherwise the later volumes go on for the
        // batches before it (one of them may fail there first), as `search_window` of
        // RPSTBLASTN does.
        let live = failed.as_ref().map_or(spans.len(), |(batch, _)| *batch);
        let mut traced = Vec::new();
        for (number, batch) in engine
            .search_volume(
                &files,
                &ProfileLetters::new(setup, volume, starts[index]),
                pool,
                &workers,
            )?
            .batches(&spans)
            .into_iter()
            .take(live)
            .enumerate()
        {
            match batch {
                Ok(results) => traced.extend(results),
                Err(error) => {
                    failed = Some((number, error));
                    break;
                }
            }
        }
        // The volume's mappings and per-volume data go before the next volume is opened (and
        // before its lists are merged): at most one volume's pages are resident.
        drop(engine);
        drop(masked);
        let aux_karlin_k = files.aux.karlin_k.clone();
        drop(files);
        if let Some((batch, _)) = &failed {
            // Nothing of this volume is merged: the run ends with the failing batch's error, and
            // the later volumes only look for an earlier batch's error.
            if *batch == 0 {
                break;
            }
            continue;
        }
        // The volume's local ordinal ids become the database's.
        let start = starts[index];
        let mut lists: Vec<Option<BlastpHitList>> = Vec::with_capacity(traced.len());
        let mut k_marks = Vec::with_capacity(traced.len());
        for traced in traced {
            let mut list = traced.hit_list;
            if let Some(list) = list.as_mut() {
                shift_ordinal_ids(list, start);
                // (The 100-slot reservation of NCBI's pointer array is 4 KB of `BlastpHspList`s.)
                list.hsplist_array.truncate(list.hsplist_count);
                list.hsplist_array.shrink_to_fit();
            }
            lists.push(list);
            k_marks.push(traced.k_mark);
        }
        if position == 0 {
            // The K of the footers: volume 0's marks and `.aux` Ks (`thread_results[0]`).
            gapped_k = batch_gapped_k(
                &run.batches,
                &k_marks,
                &aux_karlin_k,
                run.gapped_params.k,
                cbs > 0,
            );
            for (slot, list) in hit_lists.iter_mut().zip(lists) {
                *slot = list;
            }
        } else {
            // Merged into the running lists at once: no later volume's lists are held.
            merge_volume(hit_lists, lists);
        }
    }
    drop(workers);
    // The reports look the profiles up by ordinal id (`report_ordinal_ids`), batch after batch as
    // NCBI's loop reports them, up to a batch whose search failed: a profile that they would not
    // find among those they report (the first `-max_target_seqs` of each query's merged list)
    // ends the run at its batch; one past that cut leaves the list, its HSPs counted
    // (`RpsSearched::unlisted`).
    let live = failed
        .as_ref()
        .map_or(run.batches.len(), |(batch, _)| *batch);
    let mut unlisted: Vec<usize> = Vec::new();
    for (number, batch) in run.batches.iter().enumerate().take(live) {
        match report_ordinal_ids(
            db,
            &mut hit_lists[batch.clone()],
            args.max_target_seqs,
            program,
            units.searches_volumes_bare(),
        ) {
            Ok(None) => {}
            Ok(Some(dropped)) => {
                if unlisted.is_empty() {
                    unlisted = vec![0; hit_lists.len()];
                }
                unlisted[batch.clone()].copy_from_slice(&dropped.hsps);
            }
            Err(refusal) => {
                failed = Some((number, refusal));
                break;
            }
        }
    }
    if let Some((batch, error)) = failed {
        return Err(failed_batch_error(run, batch, error));
    }
    Ok(RpsSearched {
        gapped_k,
        volumes: order.len(),
        unlisted,
    })
}

/// The `EffectiveSearchSpace` option that the search of a query batch leaves in NCBI's options
/// handle, which the search strategy's next request holds, for the batch's first query `ctx`
/// (`RpsUnits::request_space` of the last unit; `-searchsp` as given, also for an invalid query,
/// whose space is 0 without it).
#[cold]
#[inline(never)]
pub(super) fn request_space(
    units: &RpsUnits,
    args: &ResolvedBlastpArgs,
    ctx: &QueryContext,
    gapped: &KarlinParams,
) -> Result<i64> {
    let program = args
        .rps_search
        .as_deref()
        .map_or("RPSBLAST", RpsSearch::program);
    match args.searchsp.unwrap_or(0) {
        0 if ctx.is_valid && ctx.aa_len > 0 => {
            units.request_space(ctx.aa_len as i64, gapped, program)
        }
        0 => Ok(0),
        searchsp => Ok(searchsp),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:193-213
/// The run's error when the query batch `batch` of `run.batches` fails (`rps_search`). NCBI reads,
/// writes the search strategy request of, searches and reports one batch before it reads the next,
/// so an error of its own in a batch's search (a profile volume's sequence file that the batch's
/// traceback reads, `ProfileLetters`: `BLAST Database error`, exit status 2) ends the run after the
/// reports of the batches before it (without the epilog) and after that batch's reader messages
/// and request. That is the error of a header file that a batch reaches (`ErrorAfterBatchReports`, named
/// by the batch's first query), whose batches before blastp's driver searches again and reports
/// alone (`search`). LOSAT's refusals end the run as they are.
#[cold]
#[inline(never)]
fn failed_batch_error(run: &RpsSearchRun<'_>, batch: usize, error: anyhow::Error) -> anyhow::Error {
    if !error.is::<crate::cli::NativeError>() {
        return error;
    }
    // The traceback fails before the batch's Seq-aligns: NCBI read none of its headers.
    crate::blastinput::app::ErrorAfterBatchReports::of_batch(run.batches[batch].start, error).into()
}

// `MAX_DBSEQ_LEN`: the longest subject (here a user query) that NCBI's engine searches whole (the
// one copy of the macro).
//
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:54
use crate::core::subject_split::MAX_DBSEQ_LEN;

#[cfg(test)]
#[path = "rps_hooks_tests.rs"]
mod tests;
