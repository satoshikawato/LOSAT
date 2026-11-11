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
use crate::algorithm::rps_search::engine::{batch_gapped_k, RpsQuery, RpsRun, RpsSearch};
use crate::algorithm::tblastx::constants::{
    GAP_TRIGGER_BIT_SCORE, X_DROP_GAPPED_FINAL, X_DROP_GAPPED_PRELIM, X_DROP_UNGAPPED_BITS,
};
use crate::algorithm::tblastx::lookup::{
    prepare_blosum62_lookup_query_for_word_size, BlastAaLookupTable, QueryContext,
};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::seqdb::rps::{RpsOpenError, RpsVolumeFiles};
use crate::stats::{BlastGumbelBlk, KarlinParams};

/// The error of the RPS files that NCBI reports itself (`BLAST engine error: <text>`, exit
/// status 3), after what the first batch's search writes before (`preliminary_search_error`).
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
/// `ADAPTIVE_CBS` changes the preliminary hit list of a search with the statistics (the stream
/// that NCBI closes early): refused, as tblastn refuses it.
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
        if std::env::var_os("ADAPTIVE_CBS").is_some() {
            bail!(
                "the environment variable ADAPTIVE_CBS, which makes NCBI BLAST+ close its preliminary stream of composition-adjusted HSPs early, is not supported by LOSAT's RPSBLAST"
            );
        }
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
    pub(super) cutoff_subject_length: i32,
    pub(super) prelim_hitlist_size: usize,
    /// The query batches (searched-query ranges) of the run.
    pub(super) batches: Vec<std::ops::Range<usize>>,
}

/// The sequence count that the statistics of an RPS search with `-dbsize` read: NCBI's
/// `x_AdjustDbSize` gives the volumes of a database the whole database's length and count, but
/// returns early where the user gave a length (or a search space), which leaves the count 0 and
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
    if args.rps_search.is_none() || args.dbsize <= 0 {
        return None;
    }
    let volumes = database.seqdb().volumes();
    if volumes.len() < 2 {
        return None;
    }
    let first = volumes
        .iter()
        .min_by(|a, b| crate::seqdb::SeqDb::compare_volume_paths(a.base(), b.base()))
        .expect("not empty");
    Some(first.num_oids() as usize)
}

/// The number of volumes that an RPS search of the database runs (`-num_threads 0` runs a thread
/// for each: NCBI's "auto"). (A volume listed twice is refused by the database layer, as for every
/// program.)
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
    let db_volumes = db.volumes();
    // The first ordinal id of each volume (the database's order), and the volumes in the order
    // of their paths as strings (a name given twice is searched once).
    let mut starts = Vec::with_capacity(db_volumes.len());
    let mut next_start = 0u32;
    for volume in db_volumes {
        starts.push(next_start);
        next_start += volume.num_oids();
    }
    let mut order: Vec<usize> = (0..db_volumes.len()).collect();
    order.sort_by(|&a, &b| {
        crate::seqdb::SeqDb::compare_volume_paths(db_volumes[a].base(), db_volumes[b].base())
    });
    let cbs = match args.comp_based_stats.mode {
        BlastpCompositionMode::NoCompositionBasedStats => 0u8,
        BlastpCompositionMode::CompositionBasedStats => 1,
        BlastpCompositionMode::CompositionMatrixAdjust => 2,
        BlastpCompositionMode::ForceFullMatrixAdjust => 3,
    };
    let gumbel = run
        .gumbel
        .expect("an RPS-BLAST search has a Gumbel block (rpsblast takes no -searchsp)");
    let lambda = run.gapped_params.lambda;
    // BLOSUM62's score block for the matrix information of the redo.
    let sbp = if cbs > 0 {
        Some(crate::algorithm::psiblast::pssm::score_block(
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:77-88
    // ```c++
    // CBlastRPSOptionsHandle::SetInitialWordOptionsDefaults()
    // {
    //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
    //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
    //     // FIXME: extend_word_method is missing
    // }
    //
    // void
    // CBlastRPSOptionsHandle::SetGappedExtensionDefaults()
    // {
    //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
    //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
    // ```
    // The X-drops in bits: `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final` or these defaults
    // (7, 15, 25); the engine converts them at each volume's scale (`RpsRun::cutoffs`,
    // `RpsRun::gapped_x_drops`).
    let x_drop_ungapped_bits = rps.x_drop_ungapped_bits.unwrap_or(X_DROP_UNGAPPED_BITS);
    let x_drop_gapped_bits = rps
        .x_drop_gapped_bits
        .unwrap_or(f64::from(X_DROP_GAPPED_PRELIM));
    let x_drop_final_bits = rps
        .x_drop_final_bits
        .unwrap_or(f64::from(X_DROP_GAPPED_FINAL));
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
    let mut first_costs: Option<(String, i32, i32)> = None;
    for (position, &index) in order.iter().enumerate() {
        let volume = &db_volumes[index];
        let num_seqs = volume.num_oids();
        let lengths = (0..num_seqs)
            .map(|oid| volume.protein_length(oid))
            .collect::<std::result::Result<Vec<u32>, _>>()
            .map_err(|error| {
                crate::algorithm::common::database::database_error(error, "RPSBLAST")
            })?;
        let to_error = |error: RpsOpenError| match error {
            RpsOpenError::Ncbi(text) => anyhow::Error::new(RpsEngineError(text)),
            RpsOpenError::Refused(text) => anyhow::anyhow!(text),
        };
        let mut files = RpsVolumeFiles::open(volume.base(), &lengths).map_err(to_error)?;
        if cbs > 0 {
            files.open_freq(volume.base()).map_err(to_error)?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:72-83
        // ```c++
        //     if (rps_info) {
        //         opts_memento->m_ScoringOpts->matrix = strdup(rps_info->aux_info.orig_score_matrix);
        //         opts_memento->m_ScoringOpts->gap_open = rps_info->aux_info.gap_open_penalty;
        //         opts_memento->m_ScoringOpts->gap_extend = rps_info->aux_info.gap_extend_penalty;
        // ```
        // The search's Karlin, Gumbel and cutoff values are those of BLOSUM62 11/1 (blastp's
        // task); a database of other costs is refused, and so are volumes whose costs differ
        // (each volume sets its own in NCBI; the statistics are one set).
        let aux = &files.aux;
        if aux.matrix != "BLOSUM62" || aux.gap_open != 11 || aux.gap_extend != 1 {
            bail!(
                "a profile database of the matrix {} with gap costs {}/{} is not supported by LOSAT's RPSBLAST",
                aux.matrix,
                aux.gap_open,
                aux.gap_extend
            );
        }
        let costs = (aux.matrix.clone(), aux.gap_open, aux.gap_extend);
        if first_costs.get_or_insert_with(|| costs.clone()) != &costs {
            bail!(
                "a profile database whose volumes differ in the matrix or the gap costs is not supported by LOSAT's RPSBLAST"
            );
        }
        // The scores of an alignment of a profile add up to at most the sum of its rows' best
        // scores; above `INT4_MAX` NCBI's 32-bit scores overflow (undefined behavior). The
        // file is read through a small buffer, not the mapping.
        let bounds = files.profile_score_bounds(volume.base()).ok_or_else(|| {
            anyhow::anyhow!(
                "a profile database whose .rps file cannot be read is not supported by LOSAT's RPSBLAST"
            )
        })?;
        if let Some(profile) = bounds.iter().position(|&bound| bound > i64::from(i32::MAX)) {
            bail!(
                "a profile database whose profile {profile} can score {} (above 2147483647, where NCBI BLAST+'s 32-bit scores overflow) is not supported by LOSAT's RPSBLAST",
                bounds[profile]
            );
        }
        // NCBI's scaled values: `(Int4)scale` multiplies the costs and cutoffs (a whole scale of
        // 1 or more here; `RpsRun::search_volume` refuses one whose products leave `Int4` and one
        // that changes Lambda).
        let scale = aux.scale;
        if !(scale >= 1.0 && scale.fract() == 0.0) {
            bail!(
                "a profile database of the scale {scale} (not a whole number of 1 or more) is not supported by LOSAT's RPSBLAST"
            );
        }
        let concat = files.rows().len();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1129-1133
        // ```c
        //     num_db_seqs = BlastSeqSrcGetNumSeqs(seq_src);
        //     dbsize = BlastSeqSrcGetTotLen(seq_src) + num_db_seqs;
        //     if (dbsize > INT4_MAX)
        //       return -3;
        // ```
        if concat as u64 + args.window_size as u64 > (i32::MAX / 2) as u64 {
            bail!("a profile volume of {concat} positions is not supported by LOSAT's RPSBLAST");
        }
        // The consensus letters at the profiles' offsets (the traceback's query letters, read
        // as stored: `Blast_TracebackGetEncoding` gives `eBlastEncodingProtein`).
        let offsets = files.offsets();
        let mut consensus = vec![0u8; concat];
        if num_seqs > 0 {
            let block = volume.read_block(0, num_seqs, 0).map_err(|error| {
                crate::algorithm::common::database::database_error(error, "RPSBLAST")
            })?;
            for oid in 0..num_seqs as usize {
                let letters = block.protein(oid);
                let start = offsets[oid] as usize;
                consensus[start..start + letters.len()].copy_from_slice(letters);
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:177-190
        // ```c
        //    case eRPSLookupTable:
        //        {
        //            BlastRPSLookupTable *lookup;
        //            Int4 alphabet_size;
        //            RPSLookupTableNew(rps_info, (BlastRPSLookupTable* *)
        //                                         (&lookup_wrap->lut));
        //
        //            /* if the alphabet size from the RPS database is too
        //               small, mask all unsupported query letters */
        //            lookup = (BlastRPSLookupTable*)(lookup_wrap->lut);
        //            alphabet_size = lookup->alphabet_size;
        //            if (alphabet_size < BLASTAA_SIZE)
        //                Blast_MaskUnsupportedAA(query, alphabet_size);
        //            break;
        //        }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1336-1346
        // ```c
        // Blast_MaskUnsupportedAA(BLAST_SequenceBlk* seq, Uint1 min_invalid)
        // {
        //     ...
        //     for (i = 0; i < length; i++) {
        //         if (sequence[i] >= min_invalid) {
        //             sequence[i] = kProtMask;
        //         }
        //     }
        // }
        // ```
        // The lookup table is made after the score block, so the queries' Karlin-Altschul values
        // keep their letters; the search (scan, extensions, traceback, rescale) reads them masked.
        const PROT_MASK: u8 = 21;
        let masked: Vec<Option<Vec<u8>>> = run
            .contexts
            .iter()
            .map(|ctx| {
                let inner = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
                (files.alphabet < crate::core::blast_util::BLASTAA_SIZE
                    && inner
                        .iter()
                        .any(|&letter| usize::from(letter) >= files.alphabet))
                .then(|| {
                    let mut copy = ctx.aa_seq.clone();
                    let end = copy.len() - 1;
                    for letter in &mut copy[1..end] {
                        if usize::from(*letter) >= files.alphabet {
                            *letter = PROT_MASK;
                        }
                    }
                    copy
                })
            })
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:932-940
        // ```c
        //     int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
        //     params->prelim_evalue = cbs_stretch*evalue;
        // ```
        let prelim_evalue = if cbs > 1 {
            5.0 * args.evalue
        } else {
            args.evalue
        };
        let engine = RpsRun {
            queries: run
                .contexts
                .iter()
                .zip(&first_query_length)
                .zip(&masked)
                .zip(&batch_kbp)
                .map(|(((ctx, &first), masked), &batch_kbp)| RpsQuery {
                    aa_seq: masked.as_deref().unwrap_or(&ctx.aa_seq),
                    valid: ctx.is_valid && ctx.aa_len > 0,
                    kbp_std: ctx.karlin_params,
                    first_query_length: first,
                    batch_kbp,
                    // The letters as read where the searched ones are masked: before SEG and the
                    // lower-case masks (`aa_seq_nomask`), else before the 26-column mask. The
                    // tabular report with a field list counts its identities on the rows itself,
                    // inside NCBI's alignment block, which the engine's count opens
                    // (`tabular::needs_rows`): it keeps that count.
                    aa_seq_nomask: if recount_identities {
                        ctx.aa_seq_nomask
                            .as_deref()
                            .or(masked.as_ref().map(|_| ctx.aa_seq.as_slice()))
                    } else {
                        None
                    },
                })
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
            x_drop_ungapped_bits,
            x_drop_gapped_bits,
            x_drop_final_bits,
            gap_trigger_bits: GAP_TRIGGER_BIT_SCORE,
            options: rps,
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
        let traced = engine.search_volume(&files, &consensus, pool)?;
        // The volume's mappings and per-volume data go before the next volume is opened (and
        // before its lists are merged): at most one volume's pages are resident.
        drop(engine);
        drop(masked);
        drop(consensus);
        let aux_karlin_k = files.aux.karlin_k.clone();
        drop(files);
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
    Ok(RpsSearched {
        gapped_k,
        volumes: order.len(),
    })
}

/// A volume's local ordinal ids as the database's: `start` is the volume's first ordinal id.
fn shift_ordinal_ids(hit_list: &mut BlastpHitList, start: u32) {
    if start == 0 {
        return;
    }
    let count = hit_list.hsplist_count;
    for list in hit_list.hsplist_array.iter_mut().take(count) {
        list.oid += start;
        for hsp in &mut list.hsps {
            hsp.s_idx += start;
        }
    }
}

/// The volumes after the first merge into the first's hit lists, query by query and volume by
/// volume, in the order of `s_MergeAlignSet`: the lists are the groups of Seq-aligns of one
/// subject (in the order of the Seq-aligns, `hit_list_to_seqalign_order`); the first Seq-align of
/// the current group of each side decides, the smaller E-value first, at equal E-values the
/// larger bit score, at equal both the final side's; a group of the input goes in whole before
/// the final group it beats, and when the final side runs out the rest of the input is
/// appended. Nothing is cut: the formatter keeps the first `max_target_seqs` subjects of the
/// merged set, and outfmt 7 counts its HSPs before the cut.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:69-137
/// ```c++
///     while(input_it != input_list.end())
///     {
///         double final_evalue;
///         double input_evalue;
///
///         (*final_it)->GetNamedScore(CSeq_align::eScore_EValue, final_evalue);
///         (*input_it)->GetNamedScore(CSeq_align::eScore_EValue, input_evalue);
///
///         if(input_evalue == final_evalue)
///         {
///             //Pulling a trick here to keep the program flow simple
///             //Replace the final evalue with input bitscore and vice versa
///             (*final_it)->GetNamedScore(CSeq_align::eScore_BitScore, input_evalue);
///             (*input_it)->GetNamedScore(CSeq_align::eScore_BitScore, final_evalue);
///         }
///
///         if(input_evalue <  final_evalue)
///         {
///             ... (the whole group of the input's subject)
///             final_list.insert(final_it, start_input_it, input_it);
///         }
///         else
///         {
///             ... (final_it past the whole group of its subject)
///             if(final_it == final_list.end())
///             {
///                 final_list.insert(final_it, input_it, input_list.end());
///                 break;
///             }
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:155-160
/// ```c++
///             if(thread_results[d]->HasAlignments())
///             {
///                 CConstRef<CSeq_align_set>  thread_align_set = thread_results[d]->GetSeqAlign();
///                 if(align_set->IsEmpty())
///                 {
///                     align_set->Set().insert(align_set->Set().begin(),
///                                         thread_align_set->Get().begin(),
///                                         thread_align_set->Get().end());
/// ```
#[cold]
#[inline(never)]
fn merge_volume(hit_lists: &mut [Option<BlastpHitList>], volume: Vec<Option<BlastpHitList>>) {
    for (slot, input) in hit_lists.iter_mut().zip(volume) {
        let Some(mut input) = input else {
            continue;
        };
        crate::algorithm::blastp::hsp::hit_list_to_seqalign_order(&mut input);
        if input.hsplist_count == 0 {
            continue;
        }
        match slot.as_mut() {
            None => {
                input.hsplist_array.shrink_to_fit();
                *slot = Some(input);
            }
            Some(final_list) => {
                crate::algorithm::blastp::hsp::hit_list_to_seqalign_order(final_list);
                if final_list.hsplist_count == 0 {
                    input.hsplist_array.shrink_to_fit();
                    *slot = Some(input);
                    continue;
                }
                merge_align_set(final_list, input);
            }
        }
    }
}

/// The E-value of a Seq-align's score list: 0 below `SMALLEST_EVALUE` (what `s_MergeAlignSet`
/// reads with `GetNamedScore`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1131 (and `:62`)
/// ```c++
/// #define SMALLEST_EVALUE 1.0e-180
/// ...
///     double evalue = (hsp->evalue < SMALLEST_EVALUE) ? 0.0 : hsp->evalue;
/// ```
fn seqalign_evalue(e_value: f64) -> f64 {
    if e_value < 1.0e-180 {
        0.0
    } else {
        e_value
    }
}

/// `s_MergeAlignSet` of one query's lists (see [`merge_volume`]).
fn merge_align_set(final_list: &mut BlastpHitList, mut input: BlastpHitList) {
    final_list.hsplist_array.truncate(final_list.hsplist_count);
    input.hsplist_array.truncate(input.hsplist_count);
    let mut incoming = input.hsplist_array.into_iter().peekable();
    let mut at = 0usize;
    while let Some(group) = incoming.peek() {
        let (final_first, input_first) = (&final_list.hsplist_array[at].hsps[0], &group.hsps[0]);
        // The Seq-align carries 0 for an E-value below `SMALLEST_EVALUE`, so those compare equal.
        let (mut final_value, mut input_value) = (
            seqalign_evalue(final_first.e_value),
            seqalign_evalue(input_first.e_value),
        );
        if input_value == final_value {
            // The trick of NCBI's loop: the bit scores take the e-values' places.
            input_value = final_first.bit_score;
            final_value = input_first.bit_score;
        }
        if input_value < final_value {
            let group = incoming.next().expect("peeked");
            final_list.hsplist_array.insert(at, group);
            at += 1;
        } else {
            at += 1;
            if at == final_list.hsplist_array.len() {
                final_list.hsplist_array.extend(incoming);
                break;
            }
        }
    }
    final_list.hsplist_count = final_list.hsplist_array.len();
    final_list.hsplist_array.shrink_to_fit();
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
