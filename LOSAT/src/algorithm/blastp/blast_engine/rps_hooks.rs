//! The stages of an rpsblast search (series PR, item PS11a) that differ from blastp's, called
//! once per run from `search_cli` and `run_resolved_in_pool` when the search is RPS-BLAST's
//! (`ResolvedBlastpArgs::rps_search`): the check of the database's frequency ratios, the queries'
//! contexts, and the whole RPS search (`algorithm::rps_search`) into the final hit lists, whose
//! subjects are the profiles' consensus sequences of the database. blastp's driver does the rest
//! (the query batches, the database's statistics, the reports). Every function is out of line
//! and cold: blastp's default path only tests `rps_search`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:129-141
//! ```c++
//!         	opts_hndl.Reset(&*m_CmdLineArgs->SetOptions(args));
//!         }
//!
//!         CheckForFreqRatioFile(m_CmdLineArgs->GetBlastDatabaseArgs()->GetDatabaseName(),
//!         		               opts_hndl, true);
//!         const CBlastOptions& opt = opts_hndl->GetOptions();
//!
//!         /*** Initialize the database ***/
//!         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
//!         CRef<CLocalDbAdapter> db_adapter;
//!         CRef<CScope> scope;
//!         InitializeSubject(db_args, opts_hndl, m_CmdLineArgs->ExecuteRemotely(),
//!                          db_adapter, scope);
//! ```

use anyhow::{bail, Result};

use super::{database, BlastpRuntimeLookup};
use crate::algorithm::blastp::args::{BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::rps_search::engine::{
    batch_gapped_k, RpsFrame, RpsQuery, RpsRun, RpsSearch, RpsWorkers,
};
use crate::algorithm::rps_search::volumes::{
    mask_unsupported_letters, merge_volume, shift_ordinal_ids, Consensus, RpsVolumeError,
    RpsVolumeSet, RpsVolumeSetup,
};
use crate::algorithm::tblastx::lookup::{
    prepare_blosum62_lookup_query_for_word_size, BlastAaLookupTable, QueryContext,
};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::stats::{BlastGumbelBlk, KarlinParams};

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:159-164
// ```c++
//     // 2. Take care of any rps information
//     if (Blast_ProgramIsRpsBlast(options->GetProgramType())) {
//         ...
//         retval->m_InternalData->m_RpsData =
//             CSetupFactory::CreateRpsStructures(rps_dbname, options);
// ```
/// The error of the RPS files that NCBI reports itself (`BLAST engine error: <text>`, exit
/// status 3): a volume's files that cannot be mapped or read, or the score block of an `.aux`
/// without a matrix name (review B of PS12a, N4). NCBI sets up every volume's files in the first
/// batch's set-up (`BlastSetupPreliminarySearchEx`, before `CheckInternalData` looks at the
/// batch's queries), so the error is the first batch's whatever its queries' validity, after
/// what that batch writes before (`preliminary_search_error`; series PR, item PS12f, review S1).
#[derive(Debug)]
pub(super) struct RpsEngineError(pub(super) String);

impl std::fmt::Display for RpsEngineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BLAST engine error: {}", self.0)
    }
}

impl std::error::Error for RpsEngineError {}

/// `CheckForFreqRatioFile`: with composition-based statistics, each database name of `-db` whose
/// first volume has no `.freq` file gets NCBI's warning and the search runs without them
/// (`args.comp_based_stats` becomes mode 0). Only the first volume of each name is looked at: a
/// later volume without `.freq` fails inside the search (`RpsVolumeFiles::open_freq`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:807-840
/// ```c++
///     bool use_cbs = (opt_handle->GetOptions().GetCompositionBasedStats() == eNoCompositionBasedStats) ? false : true;
///     if(use_cbs) {
///         vector<string> db;
///         NStr::Split(rps_dbname, " ", db);
///         list<string> failed_db;
///         for (unsigned int i=0; i < db.size(); i++) {
///     	    string path;
///     	    try {
///                 vector<string> dbpath;
///        	        CSeqDB::FindVolumePaths(db[i], CSeqDB::eProtein, dbpath);
///                 path = *dbpath.begin();
///             } catch (const CSeqDBException& e) {
///                  NCBI_RETHROW(e, CBlastException, eRpsInit,
///                               "Cannot retrieve path to RPS database");
///             }
///
///     	    CFile f(path + ".freq");
///             if(!f.Exists()) {
///             	failed_db.push_back(db[i]);
///             }
///
///         }
///         if(!failed_db.empty()) {
///         	opt_handle->SetOptions().SetCompositionBasedStats(eNoCompositionBasedStats);
///         	string all_failed = NStr::Join(failed_db, ", ");
///         	string prog_str = isRpsblast ? "RPSBLAST": "DELTABLAST";
///         	string msg = all_failed + " contain(s) no freq ratios " \
///                      	 + "needed for composition-based statistics.\n" \
///                      	 + prog_str + " will be run without composition-based statistics.";
///         	ERR_POST(Warning << msg);
///         }
///
///     }
/// ```
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
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => bail!("{reason} is not supported by LOSAT's RPSBLAST"),
    };
    let path = crate::seqdb::SearchPath::from_environment(registry).map_err(|_| {
        crate::blastinput::app::engine_error("Cannot retrieve path to RPS database")
    })?;
    let mut failed: Vec<&str> = Vec::new();
    for db in name.split(' ').filter(|db| !db.is_empty()) {
        let volumes =
            crate::seqdb::SeqDb::volume_and_alias_paths(db, crate::seqdb::SeqType::Protein, &path)
                .ok()
                .map(|(volumes, _)| volumes)
                .filter(|volumes| !volumes.is_empty())
                .ok_or_else(|| {
                    crate::blastinput::app::engine_error("Cannot retrieve path to RPS database")
                })?;
        if !std::path::Path::new(&format!("{}.freq", volumes[0])).exists() {
            failed.push(db);
        }
    }
    if failed.is_empty() {
        return Ok(());
    }
    args.comp_based_stats.mode =
        crate::algorithm::blastp::args::BlastpCompositionMode::NoCompositionBasedStats;
    writeln!(
        diagnostics,
        "Warning: [rpsblast] {} contain(s) no freq ratios needed for composition-based statistics.\nRPSBLAST will be run without composition-based statistics.",
        failed.join(", ")
    )?;
    Ok(())
}

/// What rpsblast sets in the options handle that the search strategy's requests hold
/// (`report::archive::OptionsHandle::set_search_sets`): before the batches, the composition-based
/// statistics that `check_freq_ratio_file` turned off (NCBI's `SetCompositionBasedStats` in place:
/// `CompositionBasedStats 0` in every request, NCBI BLAST+ 2.17.0 `p9.rps.err.nofreq`); at each
/// batch's set-up, the database's length and `.aux` costs (`rps_search::volumes::search_sets`).
#[cold]
#[inline(never)]
pub(super) fn note_strategy_sets(
    strategy: &crate::report::archive::strategy::StrategyExport,
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
    // volumes up again, `rps_search`).
    let volumes = RpsVolumeSet::of(database.seqdb());
    let costs = RpsVolumeSetup::new("RPSBLAST", args.window_size as i32).costs(volumes.first());
    strategy.options().set_search_sets(
        before_batches,
        crate::algorithm::rps_search::volumes::search_sets(
            database.seqdb(),
            costs,
            args.searchsp,
            Some(args.dbsize),
        ),
    )
}

/// The queries' contexts (their letters and `kbp_std`, the composition with BLOSUM62, and
/// validity) and an empty lookup table: the RPS search scans with the database's prebuilt table.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2746-2748
/// ```c
///    Boolean check_ideal =
///       (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
///        program == eBlastTypeRpsTblastn);
/// ```
#[cold]
#[inline(never)]
pub(super) fn rps_contexts(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    let (_, _, contexts) =
        prepare_blosum62_lookup_query_for_word_size(query_frames, args.word_size, false);
    let table = BlastAaLookupTable {
        backbone: Vec::new(),
        overflow: Vec::new(),
        pv: Vec::new(),
        frame_bases: Vec::new(),
        num_contexts: 0,
        query_length: 0,
        word_length: args.word_size as i32,
        alphabet_size: crate::core::blast_util::BLASTAA_SIZE as i32,
        charsize: 5,
        mask: 0,
        longest_chain: 0,
        threshold: 0,
        row_max: Vec::new(),
    };
    (BlastpRuntimeLookup::Aa(table), contexts)
}

/// The per-run values of the search that `run_resolved_in_pool` computed for blastp.
pub(super) struct RpsSearchRun<'a> {
    pub(super) contexts: &'a [QueryContext],
    /// `kbp_gap_std`: BLOSUM62 11/1's table values.
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) gumbel: Option<&'a BlastGumbelBlk>,
    /// Each query's effective search space over the whole database or `-searchsp`
    /// (`exact_search_spaces`; `None` where the database has no letters); the per-volume values
    /// of `-dbsize` with several volumes are made by `volume_search_space`.
    pub(super) search_spaces: &'a [Option<i64>],
    pub(super) cutoff_subject_length: i32,
    pub(super) prelim_hitlist_size: usize,
    /// The query batches (searched-query ranges) of the run.
    pub(super) batches: Vec<std::ops::Range<usize>>,
}

/// The sequence count that the statistics of an RPS search with `-dbsize` or `-searchsp` read:
/// NCBI's `x_AdjustDbSize` gives the volumes of a database the whole database's length and
/// count, but returns early where the user gave a length or a search space (`-searchsp` above 0;
/// the search space itself is then the user's, only the length adjustment reads the count),
/// which leaves the count 0 and
/// each volume's search reads its own: volume 0's, the first of the volumes in path order, is
/// the one whose effective search space the footers print (`thread_results[0]`'s ancillary
/// data). `None` where the database's count is the one to use.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:351-354
/// ```c++
/// 	if(m_opt_handle->GetOptions().GetEffectiveSearchSpace()!= 0)
/// 		return;
///
/// 	if(m_opt_handle->GetOptions().GetDbLength()!= 0)
/// 		return;
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-727
/// ```c
///    if (eff_len_options->db_length > 0)
///       db_length = eff_len_options->db_length;
///    else
///       db_length = eff_len_params->real_db_length;
///    ...
///    if (eff_len_options->dbseq_num > 0)
///       db_num_seqs = eff_len_options->dbseq_num;
///    else
///       db_num_seqs = eff_len_params->real_num_seqs;
/// ```
#[cold]
#[inline(never)]
pub(super) fn dbsize_volume_num_seqs(
    args: &ResolvedBlastpArgs,
    database: &database::BlastpDatabase,
) -> Option<usize> {
    if args.rps_search.is_none() || (args.dbsize <= 0 && args.searchsp.unwrap_or(0) == 0) {
        return None;
    }
    let volumes = RpsVolumeSet::of(database.seqdb());
    if volumes.order.len() < 2 {
        return None;
    }
    Some(volumes.first().num_oids() as usize)
}

/// Whether the volumes of rpsblast's search have their own effective search spaces
/// (`rps_search::volumes::per_volume_search_spaces`).
pub(super) fn per_volume_search_spaces(args: &ResolvedBlastpArgs, volumes: usize) -> bool {
    crate::algorithm::rps_search::volumes::per_volume_search_spaces(
        args.dbsize,
        args.searchsp.unwrap_or(0),
        volumes,
    )
}

/// One query's effective search space in a volume of `num_seqs` sequences with the database
/// length `db_length` (`-dbsize`): `(db_length - num_seqs * adj) * (query_length - adj)` with
/// NCBI's length adjustment, and its floor of 1 on the effective database length. A product
/// that does not fit NCBI's `Int8` is refused (NCBI's arithmetic overflows).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:826-846
/// ```c
///          if (effective_search_space == 0) {
///         	 Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
///
///         	 // Just in case effective_db_length < 0
///         	 if (effective_db_length <= 0)
///         		 effective_db_length = 1;
///
///              effective_search_space = effective_db_length *
///                              (query_length - length_adjustment);
///          }
/// ```
/// (Not `database::blastp_database_search_space`, which blastp's default path inlines.)
#[cold]
#[inline(never)]
pub(super) fn volume_search_space(
    query_length: usize,
    db_length: i64,
    num_seqs: i64,
    gapped: &KarlinParams,
) -> Result<i64> {
    let overflow = || {
        anyhow::anyhow!(
            "a database length of {db_length} letters with a query of {query_length} letters overflows NCBI BLAST+'s 64-bit effective search space, which is not supported by LOSAT's RPSBLAST"
        )
    };
    let query_length = i64::try_from(query_length).map_err(|_| overflow())?;
    let adjustment = crate::stats::length_adjustment::compute_length_adjustment_ncbi(
        query_length,
        db_length,
        num_seqs,
        gapped,
    )
    .length_adjustment;
    adjusted_search_space(query_length, db_length, num_seqs, adjustment).ok_or_else(overflow)
}

/// `(db_length - num_seqs * adjustment) * (query_length - adjustment)` in NCBI's `Int8`, with the
/// effective database length raised to 1 when it is 0 or below; `None` where a product leaves
/// `Int8`.
fn adjusted_search_space(
    query_length: i64,
    db_length: i64,
    num_seqs: i64,
    adjustment: i64,
) -> Option<i64> {
    let mut effective_db_length = num_seqs
        .checked_mul(adjustment)
        .and_then(|adjusted| db_length.checked_sub(adjusted))?;
    if effective_db_length <= 0 {
        effective_db_length = 1;
    }
    effective_db_length.checked_mul(query_length - adjustment)
}

/// LOSAT's refusal of an RPS-BLAST database whose alias files give an OID range or a key of the
/// statistics or the shortest length (`SeqDb::alias_ranges_or_statistics`): NCBI's RPS search
/// opens each volume alone and searches every profile of its lookup files
/// (rpsblast_local.cpp:341-400), which LOSAT has not traced for these keys (line RG, D1a-2).
pub(super) fn check_alias_keys(database: &database::BlastpDatabase, program: &str) -> Result<()> {
    if database.seqdb().alias_ranges_or_statistics() {
        return Err(crate::algorithm::common::database::refusal(
            "an alias file's FIRST_OID, LAST_OID, STATS_NSEQ, STATS_TOTLEN or MIN_SEQ_LENGTH over an RPS-BLAST database",
            program,
        ));
    }
    Ok(())
}

/// The number of volumes that an RPS search of the database runs (`-num_threads 0` runs a thread
/// for each: NCBI's "auto"). A volume that the alias files reach twice is one volume of the
/// database's set (`SeqDb`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:419-424
/// ```c++
///    	if((kAutoThreadedSearch == m_num_of_threads) ||
///    	  (m_num_of_threads > m_rps_databases.size()))
///    	{
///    		//Default num of thread : a thread for each db
///    		m_num_of_threads = static_cast<unsigned int>(m_rps_databases.size());
///    	}
/// ```
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
}

/// The RPS search of every query of the run against each volume of the database in turn: the
/// final hit lists of the volumes merged into `hit_lists` (by query; their subjects named by the
/// database's ordinal ids; with several volumes longer than `-max_target_seqs` and not to be
/// sorted or cut by blastp's driver before the formatter's cut), and the gapped K of each query's
/// outfmt 0 footer (volume 0's).
///
/// Every volume is searched alone, the queries in parallel, and its files are unmapped before
/// the next volume's are mapped. NCBI searches the volumes in the sorted order of their paths
/// (a `set<string>`, `FindVolumePaths`), with the statistics of the whole database
/// (`x_AdjustDbSize`: the database's length and sequence count, which `run` already holds), a
/// per-volume hit list size, and merges each query's results into volume 0's (`s_MergeAlignSet`,
/// `merge_volume`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:341-346
/// ```c++
/// 	CSeqDB::FindVolumePaths(db, CSeqDB::eProtein, m_rps_databases, NULL, true, true);
/// 	m_num_of_dbs = static_cast<unsigned int>(m_rps_databases.size());
/// 	if( 1 == m_num_of_dbs)
/// 	{
/// 		m_num_of_threads = kDisableThreadedSearch;
/// 	}
/// ```
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
    let Some(database) = database else {
        bail!("a search without -db is not supported by LOSAT's RPSBLAST");
    };
    let db = database.seqdb();
    let RpsVolumeSet {
        volumes: db_volumes,
        starts,
        order,
    } = RpsVolumeSet::of(db);
    let cbs = match args.comp_based_stats.mode {
        BlastpCompositionMode::NoCompositionBasedStats => 0u8,
        BlastpCompositionMode::CompositionBasedStats => 1,
        BlastpCompositionMode::CompositionMatrixAdjust => 2,
        BlastpCompositionMode::ForceFullMatrixAdjust => 3,
    };
    // The Gumbel block, none with `-searchsp` (at any value) or `OLD_FSC`.
    let gumbel = run.gumbel;
    let per_volume_spaces = per_volume_search_spaces(args, db_volumes.len());
    let lambda = run.gapped_params.lambda;
    // BLOSUM62's score block for the matrix information of the redo.
    let sbp = if cbs > 0 {
        Some(crate::algorithm::pssm_search::in_pssm::score_block(
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
    // NCBI's reports count an HSP's identities on the sequences of their scope (outfmt 0 always,
    // the tabular formats inside their alignment block), not on the masked query of the
    // traceback's count (`rps_search::traceback::formatter_identities`). LOSAT's outfmt 0 and
    // default tabular writers print the HSP's count, so the engine counts again for them; the
    // tabular writer of a field list or outfmt 10 counts on the rows itself and reads the
    // traceback's count for the block's gate (`tabular::needs_rows`).
    let recount_identities = super::super::tabular::custom_tabular_spec(&args.outfmt).is_none();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250
    // ```c
    //     if (backup->offset + MAX_DBSEQ_LEN <
    //         backup->hard_ranges[backup->hm_index].right) {
    //
    //         subject->length = MAX_DBSEQ_LEN;
    //         backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
    // ```
    // A query that NCBI's engine would cut into chunks (`s_GetNextSubjectChunk`) is refused.
    if let Some(long) = run.contexts.iter().find(|ctx| ctx.aa_len > MAX_DBSEQ_LEN) {
        bail!(
            "a query of {} letters (NCBI BLAST+ searches a query above {MAX_DBSEQ_LEN} letters in chunks) is not supported by LOSAT's RPSBLAST",
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
    let mut setup = RpsVolumeSetup::new("RPSBLAST", args.window_size as i32);
    // The workers' scratch, kept from volume to volume (`RpsWorkers`).
    let workers = RpsWorkers::new(pool.threads());
    // The first failing batch (its number and error) of the volumes searched so far.
    let spans: Vec<usize> = run.batches.iter().map(|batch| batch.len()).collect();
    let mut failed: Option<(usize, anyhow::Error)> = None;
    for (position, &index) in order.iter().enumerate() {
        let volume = &db_volumes[index];
        let num_seqs = volume.num_oids();
        let files = match setup.open(volume, cbs > 0)? {
            Ok(files) => files,
            Err(RpsVolumeError::Files(text)) => return Err(RpsEngineError(text).into()),
            // N2 of review A (PS12a): the score block of each query of the first batch fails.
            Err(RpsVolumeError::NoMatrix) => {
                let queries = run.batches.first().map_or(0, |batch| batch.len());
                return Err(RpsEngineError("Error: Unknown error code -1 ".repeat(queries)).into());
            }
        };
        let masked: Vec<Option<Vec<u8>>> = run
            .contexts
            .iter()
            .map(|ctx| mask_unsupported_letters(&ctx.aa_seq, files.alphabet))
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:932-940
        // ```c
        //     int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
        //     params->prelim_evalue = cbs_stretch*evalue;
        // ```
        // (The stretch is NCBI's only inside the Gumbel block's branch: blast_parameters.c:922-945.)
        let prelim_evalue = if cbs > 1 && gumbel.is_some() {
            5.0 * args.evalue
        } else {
            args.evalue
        };
        // Each query's effective search space in this volume.
        let search_spaces = run
            .contexts
            .iter()
            .zip(run.search_spaces)
            .map(|(ctx, &whole)| {
                if !(ctx.is_valid && ctx.aa_len > 0) {
                    Ok(0)
                } else if per_volume_spaces {
                    volume_search_space(
                        ctx.aa_len,
                        args.dbsize,
                        i64::from(num_seqs),
                        run.gapped_params,
                    )
                } else {
                    Ok(whole.unwrap_or(0))
                }
            })
            .collect::<Result<Vec<i64>>>()?;
        let engine = RpsRun {
            queries: run
                .contexts
                .iter()
                .zip(&first_query_length)
                .zip(&masked)
                .zip(&batch_kbp)
                .zip(&search_spaces)
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
            gumbel,
            cutoff_subject_length: run.cutoff_subject_length,
            evalue: args.evalue,
            prelim_evalue,
            cbs,
            sbp: sbp.as_ref(),
            hitlist_size: args.traceback_hitlist_size(),
            prelim_hitlist_size: run.prelim_hitlist_size,
            window: args.window_size as i32,
            x_drop_ungapped_bits: x_drops.ungapped,
            x_drop_gapped_bits: x_drops.gapped,
            x_drop_final_bits: x_drops.final_,
            gap_trigger_bits: x_drops.gap_trigger,
            options: rps,
            program: "RPSBLAST",
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873
        // ```c
        //     if (hit_params->link_hsp_params) {
        //         status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
        //                   subject->length, gap_align->sbp, hit_params->link_hsp_params,
        //                   score_options->gapped_calculation);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:559-561
        // ```c
        //       query_context = hp_start->next->hsp->context;
        //       length_adjustment = query_info->contexts[query_context].length_adjustment;
        //       query_length = query_info->contexts[query_context].query_length;
        // ```
        // With sum statistics (`-sum_stats true`) NCBI's RPS search links the HSPs of each
        // profile with the user query's one-context `query_info` while the HSPs carry the
        // profile's context (0 to the number of profiles), so `BLAST_LinkHsps` reads past its
        // arrays (valgrind: invalid reads; a segmentation fault on some queries; e-values of 0.0
        // otherwise): undefined behaviour, refused where the first volume's search would start,
        // after NCBI's errors of the options, the database and its files.
        if rps.sum_stats {
            bail!(
                "-sum_stats true (NCBI BLAST+ 2.17.0's rpsblast links HSPs with the profiles' contexts as query indexes and reads past its arrays: its output depends on memory contents or it crashes) is not supported by LOSAT's RPSBLAST"
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
                &Consensus::new(&setup, volume, &files),
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
    if let Some((batch, error)) = failed {
        return Err(failed_batch_error(run, batch, error));
    }
    Ok(RpsSearched {
        gapped_k,
        volumes: order.len(),
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:193-213
// ```c++
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//
//             SaveSearchStrategy(args, m_CmdLineArgs, queries, opts_hndl);
//             ...
//             	CLocalRPSBlast  local_search (query_batch, db_args->GetDatabaseName(), opts_hndl, args[kArgNumThreads].AsInteger() );
//             	results = local_search.Run();
// ```
/// The run's error when the query batch `batch` of `run.batches` fails (`rps_search`). NCBI reads,
/// writes the search strategy request of, searches and reports one batch before it reads the next,
/// so an error of its own in a batch's search (a profile volume's sequence file that the batch's
/// traceback reads, `Consensus`: `BLAST Database error`, exit status 2) ends the run after the
/// reports of the batches before it (without the epilog) and after that batch's reader messages
/// and request. That is the error of a header file that a batch reaches (`HeaderReadError`, named
/// by the batch's first query), whose batches before blastp's driver searches again and reports
/// alone (`search`). LOSAT's refusals end the run as they are.
#[cold]
#[inline(never)]
fn failed_batch_error(run: &RpsSearchRun<'_>, batch: usize, error: anyhow::Error) -> anyhow::Error {
    if !error.is::<crate::cli::NativeError>() {
        return error;
    }
    database::HeaderReadError {
        query: run.batches[batch].start,
        error,
    }
    .into()
}

/// `MAX_DBSEQ_LEN`: the longest subject (here a user query) that NCBI's engine searches whole.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:54
/// ```c
/// #define MAX_DBSEQ_LEN 5000000
/// ```
const MAX_DBSEQ_LEN: usize = 5_000_000;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::hsp::{BlastpHsp, BlastpHspList};

    /// A list of one subject whose first Seq-align has the given E-value and bit score.
    fn list(oid: u32, e_value: f64, bit_score: f64) -> BlastpHspList {
        let mut hsp = BlastpHsp::from_hit(crate::common::Hit {
            identity: 0.0,
            length: 1,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 1,
            s_start: 1,
            s_end: 1,
            e_value,
            bit_score,
            num_ident: 0,
            query_frame: 0,
            query_length: 0,
            q_idx: 0,
            s_idx: oid,
            raw_score: 10,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 0,
        });
        hsp.e_value = e_value;
        hsp.bit_score = bit_score;
        BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![hsp],
            best_evalue: e_value,
        }
    }

    fn hit_list(lists: Vec<BlastpHspList>) -> BlastpHitList {
        let mut hit_list = BlastpHitList::new(100);
        hit_list.hsplist_count = lists.len();
        hit_list.hsplist_array = lists;
        hit_list
    }

    fn oids(hit_list: &BlastpHitList) -> Vec<u32> {
        hit_list.hsplist_array[..hit_list.hsplist_count]
            .iter()
            .map(|list| list.oid)
            .collect()
    }

    fn blosum62_gapped() -> KarlinParams {
        KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        }
    }

    // blast_setup.c:826-846 (quoted at `volume_search_space`): each volume's space reads the
    // given length and its own sequence count; an effective database length of 0 or below is 1.
    #[test]
    fn a_volumes_search_space_reads_its_own_sequence_count() {
        use crate::stats::length_adjustment::compute_length_adjustment_ncbi;
        let gapped = blosum62_gapped();
        for (query, db_length, seqs) in [
            (302usize, 5_000i64, 1i64),
            (302, 5_000, 3),
            (40, 100_000, 7),
        ] {
            let adjustment = compute_length_adjustment_ncbi(query as i64, db_length, seqs, &gapped)
                .length_adjustment;
            assert_eq!(
                volume_search_space(query, db_length, seqs, &gapped).unwrap(),
                (db_length - seqs * adjustment) * (query as i64 - adjustment),
                "{query} {db_length} {seqs}"
            );
        }
        let one = volume_search_space(302, 5_000, 1, &gapped).unwrap();
        let three = volume_search_space(302, 5_000, 3, &gapped).unwrap();
        assert_ne!(one, three);
        // The floor: an effective database length of 0 or below is 1 (NCBI's length adjustment
        // keeps it above 0, so the arithmetic is tried on given adjustments).
        assert_eq!(adjusted_search_space(302, 100, 10, 10), Some(292));
        assert_eq!(adjusted_search_space(302, 100, 10, 11), Some(291));
        assert_eq!(adjusted_search_space(302, 100, 10, 9), Some(10 * 293));
        // A product past `Int8` is refused.
        assert!(volume_search_space(302, i64::MAX, 1, &gapped)
            .unwrap_err()
            .to_string()
            .contains("overflows NCBI BLAST+'s 64-bit effective search space"));
    }

    #[test]
    fn volumes_merge_by_the_first_alignment_of_each_subject() {
        // Final: 1 (1e-5), 2 (1e-3), 3 (1e-1). Input: 10 (1e-6), 11 (1e-3, a larger bit score
        // than 2), 12 (1e-3, the same bit score as 2 or less), 13 (5).
        let mut lists = vec![Some(hit_list(vec![
            list(1, 1e-5, 50.0),
            list(2, 1e-3, 40.0),
            list(3, 1e-1, 30.0),
        ]))];
        let input = vec![Some(hit_list(vec![
            list(10, 1e-6, 55.0),
            list(11, 1e-3, 41.0),
            list(12, 1e-3, 40.0),
            list(13, 5.0, 20.0),
        ]))];
        merge_volume(&mut lists, input);
        // 10 before 1; at equal E-values 11 (larger bit score) before 2, but 12 (equal) after 2;
        // 13 goes after 3 (the final side runs out: the rest is appended).
        assert_eq!(
            oids(lists[0].as_ref().unwrap()),
            vec![10, 1, 11, 2, 12, 3, 13]
        );
    }

    #[test]
    fn a_volume_without_hits_changes_nothing_and_an_empty_first_volume_takes_the_next() {
        let mut lists = vec![Some(hit_list(vec![list(1, 1e-5, 50.0)])), None];
        let input = vec![
            None,
            Some(hit_list(vec![list(10, 1e-6, 55.0), list(11, 1e-2, 9.0)])),
        ];
        merge_volume(&mut lists, input);
        assert_eq!(oids(lists[0].as_ref().unwrap()), vec![1]);
        assert_eq!(oids(lists[1].as_ref().unwrap()), vec![10, 11]);
    }

    #[test]
    fn the_input_is_cut_nowhere_and_volumes_merge_in_order() {
        let mut lists = vec![Some(hit_list(vec![list(1, 1e-2, 10.0)]))];
        let second = vec![Some(hit_list(vec![
            list(10, 1e-3, 10.0),
            list(11, 1e-1, 10.0),
        ]))];
        let third = vec![Some(hit_list(vec![list(20, 1e-2, 10.0)]))];
        merge_volume(&mut lists, second);
        merge_volume(&mut lists, third);
        // 20 ties with 1 on both values: the final side's (already merged) group goes first.
        assert_eq!(oids(lists[0].as_ref().unwrap()), vec![10, 1, 20, 11]);
    }

    #[test]
    fn volumes_are_ordered_by_file_name_then_path() {
        use crate::seqdb::SeqDb;
        use std::path::Path;
        let mut paths = ["zd/x", "ad/y", "prof_rps2v.01", "b/x", "a/x"]
            .map(Path::new)
            .to_vec();
        paths.sort_by(|a, b| SeqDb::compare_volume_paths(a, b));
        let sorted: Vec<&str> = paths.iter().map(|p| p.to_str().unwrap()).collect();
        // prof_rps2v.01, then x (a/x, b/x, zd/x by path), then y.
        assert_eq!(sorted, ["prof_rps2v.01", "a/x", "b/x", "zd/x", "ad/y"]);
    }

    #[test]
    fn e_values_below_the_smallest_compare_equal_and_the_bit_score_decides() {
        // 1e-189 (535 bits) and 1e-190 (534 bits) are both 0 in the Seq-aligns: the larger bit
        // score first, whichever has the smaller E-value.
        let mut lists = vec![Some(hit_list(vec![
            list(1, 1e-189, 535.0),
            list(2, 1e-100, 300.0),
        ]))];
        let input = vec![Some(hit_list(vec![
            list(10, 1e-190, 534.0),
            list(11, 1e-190, 536.0),
        ]))];
        merge_volume(&mut lists, input);
        // 10 (534) after 1 (535); 11 (536) is compared with 2 (1e-100): before it.
        assert_eq!(oids(lists[0].as_ref().unwrap()), vec![1, 10, 11, 2]);
        let mut lists = vec![Some(hit_list(vec![list(1, 1e-189, 535.0)]))];
        let input = vec![Some(hit_list(vec![list(10, 1e-190, 536.0)]))];
        merge_volume(&mut lists, input);
        assert_eq!(oids(lists[0].as_ref().unwrap()), vec![10, 1]);
    }
}
