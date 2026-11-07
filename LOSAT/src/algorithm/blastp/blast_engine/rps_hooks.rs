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
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
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
/// first volume has no `.freq` file gets NCBI's warning and the search runs without them; a
/// database with frequency ratios asks for the statistics that LOSAT runs in item PS11b (refused
/// until then).
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
    args: &ResolvedBlastpArgs,
    rps: &RpsSearch,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    diagnostics: &mut dyn std::io::Write,
) -> Result<()> {
    if rps.requested_cbs == 0 {
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
        bail!(
            "-comp_based_stats {} (composition-based statistics with the database's frequency ratios) is not supported by LOSAT's RPSBLAST yet",
            rps.requested_cbs
        );
    }
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

/// The RPS search of every query of the run against the database's one volume: the final hit
/// lists into `hit_lists` (by query; their subjects named by ordinal id, the profiles' local
/// ids of the one volume), and the gapped K of each query's outfmt 0 footer.
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
) -> Result<Vec<f64>> {
    let Some(database) = database else {
        bail!("a search without -db is not supported by LOSAT's RPSBLAST");
    };
    let db = database.seqdb();
    let volumes = db.volumes();
    if volumes.len() != 1 {
        bail!(
            "a profile database of {} volumes is not supported by LOSAT's RPSBLAST yet",
            volumes.len()
        );
    }
    // The one volume's local ordinal ids are the database's.
    let volume = &volumes[0];
    let num_seqs = volume.num_oids();
    let lengths = (0..num_seqs)
        .map(|oid| volume.protein_length(oid))
        .collect::<std::result::Result<Vec<u32>, _>>()
        .map_err(|error| crate::algorithm::common::database::database_error(error, "RPSBLAST"))?;
    let files = RpsVolumeFiles::open(volume.base(), &lengths).map_err(|error| match error {
        RpsOpenError::Ncbi(text) => anyhow::Error::new(RpsEngineError(text)),
        RpsOpenError::Refused(text) => anyhow::anyhow!(text),
    })?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:72-83
    // ```c++
    //     if (rps_info) {
    //         opts_memento->m_ScoringOpts->matrix = strdup(rps_info->aux_info.orig_score_matrix);
    //         opts_memento->m_ScoringOpts->gap_open = rps_info->aux_info.gap_open_penalty;
    //         opts_memento->m_ScoringOpts->gap_extend = rps_info->aux_info.gap_extend_penalty;
    // ```
    // The search's Karlin, Gumbel and cutoff values are those of BLOSUM62 11/1 (blastp's task);
    // a database of other costs is refused.
    let aux = &files.aux;
    if aux.matrix != "BLOSUM62" || aux.gap_open != 11 || aux.gap_extend != 1 {
        bail!(
            "a profile database of the matrix {} with gap costs {}/{} is not supported by LOSAT's RPSBLAST",
            aux.matrix,
            aux.gap_open,
            aux.gap_extend
        );
    }
    // NCBI's scaled values: `(Int4)scale` multiplies the costs and cutoffs (a whole scale of 1
    // or more here; `RpsRun::search_volume` refuses one whose products leave `Int4` and one that
    // changes Lambda).
    let scale = aux.scale;
    let lambda = run.gapped_params.lambda;
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
    let gumbel = run
        .gumbel
        .expect("an RPS-BLAST search has a Gumbel block (rpsblast takes no -searchsp)");
    // The consensus letters at the profiles' offsets (the traceback's query letters, read as
    // stored: `Blast_TracebackGetEncoding` gives `eBlastEncodingProtein`).
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
    let mut first_query_length = vec![0usize; run.contexts.len()];
    for batch in &run.batches {
        if let Some(first) = run.contexts.get(batch.start) {
            for length in &mut first_query_length[batch.clone()] {
                *length = first.aa_len;
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-463
    // ```c
    //    if (sbp->kbp_gap) {
    //       double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
    //       params->gap_x_dropoff = (Int4)
    //           (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
    //       ...
    //       params->gap_x_dropoff_final = (Int4)
    //           MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda, params->gap_x_dropoff);
    //    }
    // ```
    let x_drop_gapped = (X_DROP_GAPPED_PRELIM as f64 * std::f64::consts::LN_2 / lambda) as i32;
    let x_drop_final =
        ((X_DROP_GAPPED_FINAL as f64 * std::f64::consts::LN_2 / lambda) as i32).max(x_drop_gapped);
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
    let engine = RpsRun {
        queries: run
            .contexts
            .iter()
            .zip(&first_query_length)
            .zip(&masked)
            .map(|((ctx, &first), masked)| RpsQuery {
                aa_seq: masked.as_deref().unwrap_or(&ctx.aa_seq),
                valid: ctx.is_valid && ctx.aa_len > 0,
                kbp_std: ctx.karlin_params,
                first_query_length: first,
            })
            .collect(),
        gapped: *run.gapped_params,
        gumbel,
        cutoff_subject_length: run.cutoff_subject_length,
        evalue: args.evalue,
        hitlist_size: args.traceback_hitlist_size(),
        prelim_hitlist_size: run.prelim_hitlist_size,
        window: args.window_size as i32,
        x_drop_gapped,
        x_drop_final,
        x_drop_ungapped_bits: X_DROP_UNGAPPED_BITS,
        gap_trigger_bits: GAP_TRIGGER_BIT_SCORE,
        options: rps,
    };
    let traced = engine.search_volume(&files, &consensus, pool)?;
    let mut k_marks = Vec::with_capacity(traced.len());
    for (slot, traced) in hit_lists.iter_mut().zip(traced) {
        *slot = traced.hit_list;
        k_marks.push(traced.k_mark);
    }
    Ok(batch_gapped_k(
        &run.batches,
        &k_marks,
        &files.aux.karlin_k,
        run.gapped_params.k,
    ))
}

/// `MAX_DBSEQ_LEN`: the longest subject (here a user query) that NCBI's engine searches whole.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:54
/// ```c
/// #define MAX_DBSEQ_LEN 5000000
/// ```
const MAX_DBSEQ_LEN: usize = 5_000_000;
