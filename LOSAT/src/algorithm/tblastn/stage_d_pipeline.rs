//! Internal local-subject TBLASTN Stage C to D composition pipeline.
//! The public CLI remains gated while Stage D and Stage E are incomplete.

use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::wording::ENGINE_ERROR;
use std::collections::HashMap;
#[cfg(feature = "parallel")]
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{ensure, Context, Result};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::kappa::{
    convert_distinct_alignments, normalize_postredo_scores, postredo_converted_stats_with_matrix,
    postredo_engine_positives, reap_contained_postredo_hsps, ConvertedKappaHsp,
};
use super::kappa_heap::{compo_early_termination, CompoHeap, CompoHeapRecord};
use super::match_redo::{
    match_outcome, sw_batch_failed, sw_pssm_refusal, MatchRedo, SmithWatermanRedo, XDropRedo,
};
use super::search_gapped::{
    full_translation_traceback_with_matrix_and_events_with_mask_mode_owned, GappedHsp,
    PreliminaryProfile, TargetTranslation, TracebackOwnedHsp,
};
use super::stage_d_kappa_params::{
    local_extension_final_xdrop, local_kappa_redo_params, use_matrix_workspace,
};
use super::stage_d_linking::{
    compare_preliminary_lists_for_kappa, link_preliminary_hsps, link_preliminary_hsps_on,
    reap_by_evalue, LinkedHsp, LinkedHspList,
};
use super::stage_d_results::{
    merge_query_chunk, KappaHspPayload, KappaResultHitList, KappaResultList, PrelimCulling,
};
use super::stage_d_stats::{LocalParameterCall, LocalParameterOptions, LocalSubjectParameters};
use super::x_subject_store::{SubjectLetters, XPrelimSubject, XSearchSubject};
use crate::algorithm::blastn::query_split::SplitSizes;
use crate::algorithm::blastp::encoding::{
    encode_protein_query_frame_with_seg, encode_protein_query_frame_with_seg_timed,
};
use crate::algorithm::blastp::protein_query_split::{
    protein_chunk_part_in_record, protein_chunk_would_be_split, split_protein_batch,
    ProteinQueryChunk,
};
use crate::algorithm::common::database::{FetchExcluded, FetchMarks};
use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::core::composition_adjustment::adjust_scores::{
    compo_matrix_data, read_aa_composition, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    build_query_word_hashes, BlastCompoAdjustMode, BlastCompoQueryInfo, BlastCompoSequenceData,
    BlastRedoAlignParams, EMatrixAdjustRule,
};
use crate::core::composition_adjustment::smith_waterman_match::SwBatchStatus;
use crate::stats::karlin::evalue_from_raw_score;
use crate::stats::spouge::{blast_spouge_stoe, lookup_protein_gumbel_params, BlastGumbelBlk};
use crate::stats::tables::{
    lookup_protein_params, lookup_protein_params_ungapped, protein_scoring_supported, KarlinParams,
};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::aa_char_to_ncbistdaa;
use crate::utils::seg::SegParams;
use crate::utils::threading::{with_search_pool, SearchPool};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:673-675;
// c++/include/algo/blast/core/blast_options.h:128,144,163:
// gap_x_dropoff = BLAST_GAP_X_DROPOFF_PROT; /* 15 */
// gap_x_dropoff_final = BLAST_GAP_X_DROPOFF_FINAL_PROT; /* 25 */
#[derive(Clone, Copy)]
pub(super) struct LocalStageDProfile<'a> {
    pub seg: Option<&'a SegParams>,
    pub soft_masking: bool,
    pub mask_lowercase: bool,
    pub genetic_code: u8,
    pub expect_value: f64,
    pub max_target_seqs: usize,
    /// Each query's whole record and where its searched letters start in it (`-query_loc`;
    /// empty: each query is its record, from 0): the lower-case runs of the record mask a split
    /// query's chunks (`protein_chunk_part_in_record`).
    pub query_records: &'a [(&'a [u8], usize)],
    /// The numbers of the subject set that the statistics read (a database's, `-dbsize`).
    pub subject_set: SubjectSetNumbers,
    /// The Seq-aligns are requested (outfmt 8, 12; `ReportNeeds::seqalign`): the HSPs of a
    /// composition-adjusted match get NCBI's engine positives (`postredo_engine_positives`).
    pub engine_positives: bool,
    /// The `-culling_limit` value (0: no culling). A positive value makes the preliminary
    /// stage keep the HSPs of each query in NCBI's culling writer (`culling.rs`)
    /// instead of the plain hit-list collector; the pipe after the traceback is run by the
    /// finish of the search (`tblastn_app.rs`).
    pub culling_limit: i32,
    /// `-subject_besthit`: the preliminary stage filters each subject frame's list after each
    /// merge of a chunk (`search_gapped`), the finish of the search each final hit list.
    pub subject_besthit: bool,
    /// The `-max_intron_length` value (0: NCBI's default). With sum statistics it gives the
    /// largest subject gap of the linking, or turns the linking off
    /// (`core::link_hsps::translated_longest_intron`); without them nothing reads it.
    pub max_intron_length: i32,
    /// `-in_pssm` (PSI-TBLASTN): the query's PSSM, the source of the scores of the batch's one
    /// query (`score::PssmScore`) and of its Karlin-Altschul blocks (`query_set_setup`); `None`
    /// for a FASTA query batch.
    pub pssm: Option<&'a crate::algorithm::pssm_search::PsiSearch>,
    /// NCBI's fetches with `check_oid_exclusion` of a database search (`database::fetch_excluded`;
    /// `None` without a database or where they exclude nothing): the subjects they exclude leave
    /// the preliminary hit lists before the redo or the traceback.
    pub fetch_excluded: Option<&'a crate::algorithm::common::database::FetchExcluded<'a>>,
}

/// The numbers of a database (`-db`) as NCBI's setup reads them from its sequence source
/// (`SeqDb::statistics`): the total length and the number of sequences of the setup
/// (`DbStatistics::setup_total`, `setup_count`: NCBI's `Int8`, not negative, which an alias
/// file's `LENGTH` declares, so 2^32 or more also on a 32-bit target; NCBI's `Int4`, negative for
/// a negative `NSEQ`)
/// and the Gumbel block's length in nucleotides (`gumbel_db_length`: the total before the
/// fallback to OID 0, so 0 or negative for an alias file's `LENGTH` of 0 or of 2^63 or more).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct DatabaseNumbers {
    pub total_length: i64,
    pub num_seqs: i32,
    pub gumbel_db_length: i64,
    /// `BlastSeqSrcGetTotLen` (`DbStatistics::engine_db_length`): only a positive one sorts the
    /// results at the end of the traceback (`SubjectSetNumbers::sorts_results`).
    pub engine_db_length: i64,
    /// `CSeqDB::GetMinLength` in nucleotides (`SeqDb::min_length`): the alias files'
    /// `MIN_SEQ_LENGTH` or `BLAST_SEQSRC_MINLENGTH` without a list, the scanned shortest length
    /// with one (line DI-8).
    pub min_length: i32,
    /// The number of sequences of the preliminary stage's option spaces
    /// (`DbStatistics::request_count`).
    pub request_count: i32,
}

/// The numbers of the subject set behind the statistics: a database's (`database`), or,
/// without one, those of the FASTA subject set (`-subject`: the total and the number of its
/// records, its shortest record raised to `BLAST_SEQSRC_MINLENGTH`); and `-dbsize`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct SubjectSetNumbers {
    pub database: Option<DatabaseNumbers>,
    pub dbsize: Option<i64>,
    /// `-searchsp`: the effective search space when above 0.
    pub searchsp: Option<i64>,
    /// No Gumbel block (`gumbel_block`): `-searchsp` at any value, or the environment variable
    /// `OLD_FSC` (`app::old_fsc_environment`).
    pub old_fsc: bool,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:131-136
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:908-922
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-735
impl SubjectSetNumbers {
    /// `BLAST_SEQSRC_MINLENGTH` (blast_seqsrc.h): the shortest length a sequence source
    /// answers.
    const MIN_LENGTH: usize = 10;

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181
    /// The total length in nucleotides, NCBI's `Int8`: the database's, or the sum of the
    /// subjects' (each held in memory, so at most `isize::MAX`).
    pub(super) fn total_nt(&self, subjects: &SubjectLetters<'_>) -> i64 {
        match self.database {
            Some(database) => database.total_length,
            None => subjects.lengths().map(|length| length as i64).sum(),
        }
    }

    /// Whether NCBI's traceback sorts the results by e-value at its end
    /// (`DbStatistics::sorts_results`): a database's `BlastSeqSrcGetTotLen` is positive; a FASTA
    /// set that finds hits has letters.
    pub(super) fn sorts_results(&self) -> bool {
        self.database
            .is_none_or(|database| database.engine_db_length > 0)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407, 1434 (cited at
    // `seqdb::DbStatistics`)
    /// Whether NCBI's preliminary search computes the statistics of each subject on its own: a
    /// database whose `BlastSeqSrcGetTotLen` is 0 (`DbStatistics::per_subject`); a FASTA set's
    /// total is the sum of its records.
    pub(super) fn per_subject(&self) -> bool {
        self.database
            .is_some_and(|database| database.engine_db_length == 0)
    }

    /// `-dbsize` in nucleotides when above 0, which the effective lengths read in place of the
    /// database's (or a subject's) length (blast_setup.c:716-721), NCBI's `Int8`.
    pub(super) fn dbsize_nt(&self) -> Option<i64> {
        self.dbsize.filter(|&dbsize| dbsize > 0)
    }

    /// The number of sequences, NCBI's `Int4`.
    pub(super) fn num_seqs(&self, subjects: &SubjectLetters<'_>) -> i32 {
        match self.database {
            Some(database) => database.num_seqs,
            None => i32::try_from(subjects.count()).unwrap_or(i32::MAX),
        }
    }

    /// The database length of the effective lengths in nucleotides (divided by 3 there): a
    /// positive `-dbsize`, else the total `total_nt` (blast_setup.c:716-735).
    pub(super) fn eff_length_nt(&self, total_nt: i64) -> i64 {
        self.dbsize_nt().unwrap_or(total_nt)
    }

    /// The Gumbel block's database length in residues: a `-dbsize` other than 0, else a
    /// database's length before the fallback to OID 0 (`DatabaseNumbers::gumbel_db_length`) or
    /// the total `total_nt` of a FASTA set, divided by 3 (blast_setup.c:916-922; C's `Int8`
    /// division, toward 0).
    pub(super) fn gumbel_db_length(&self, total_nt: i64) -> i64 {
        match (self.dbsize, self.database) {
            (Some(dbsize), _) if dbsize != 0 => dbsize / 3,
            (_, Some(database)) => database.gumbel_db_length / 3,
            _ => total_nt / 3,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:921-923
    /// The Gumbel block of the database length `total_nt`, none with `-searchsp` or `OLD_FSC`
    /// (`old_fsc`).
    pub(super) fn gumbel_block(
        &self,
        spec: &ProteinScoringSpec,
        total_nt: i64,
    ) -> Option<BlastGumbelBlk> {
        if self.old_fsc {
            return None;
        }
        lookup_protein_gumbel_params(spec, self.gumbel_db_length(total_nt))
    }

    /// The shortest subject length in nucleotides (`BlastSeqSrcGetMinSeqLen`, NCBI's `int`), which
    /// the Gumbel block's setup divides by 3 (`stats::setup::setup_subject_length`): a database
    /// answers `CSeqDB::GetMinLength` (`SeqDb::min_length`), a FASTA set its shortest record raised
    /// to `BLAST_SEQSRC_MINLENGTH`.
    pub(super) fn shortest_subject_length(&self, subjects: &SubjectLetters<'_>) -> i32 {
        match self.database {
            Some(database) => database.min_length,
            None => subjects
                .lengths()
                .min()
                .unwrap_or(i32::MAX as usize)
                .max(Self::MIN_LENGTH) as i32,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:903-936;
// core/blast_parameters.c:302-383,902-999:
// scoring matrix, gap costs, word size, threshold and two-hit window are
// independent options supplied before the preliminary and Stage D calls.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:457-463:
// gap_x_dropoff and gap_x_dropoff_final are converted from bits using lambda.
#[derive(Clone, Copy)]
pub(super) struct LocalStageDScoring {
    pub matrix: ScoringMatrix,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub word_size: usize,
    /// The word threshold as the options hold it (`-threshold`); the lookup table reads it in its
    /// own units (`lookup_threshold`).
    pub threshold: f64,
    pub window: i32,
    pub gap_xdrop_bits: f64,
    pub final_xdrop_bits: f64,
    /// The ungapped X-drop in bits (`word_options->x_dropoff`): `-xdrop_ungap` or 7.
    pub word_xdrop_bits: f64,
    /// `-use_sw_tback`: NCBI's traceback option set to Smith-Waterman (`eSmithWatermanTbck`), the
    /// Smith-Waterman redo of each match (`match_redo::SmithWatermanRedo`).
    pub smith_waterman: bool,
}

impl Default for LocalStageDScoring {
    fn default() -> Self {
        Self {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
            word_size: 3,
            threshold: 13.0,
            window: 40,
            gap_xdrop_bits: 15.0,
            final_xdrop_bits: 25.0,
            word_xdrop_bits: 7.0,
            smith_waterman: false,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905;
// core/blast_kappa.c:3525-3565; core/blast_traceback.c:717-719:
// comparison-only boundary snapshots expose the actual Stage C HSP list,
// Kappa incoming order and ordinary posttraceback stat_length to tests.
// They are never used to produce a search result.
#[derive(Default)]
struct StageDBoundaryTrace {
    parameters: Option<LocalSubjectParameters>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:72-115:
    ungapped_karlin: Vec<KarlinParams>,
    query_validity: Vec<bool>,
    subject_nt_lengths: Vec<usize>,
    preliminary: Vec<(usize, Vec<(usize, GappedHsp)>)>,
    redo: Vec<(usize, usize, Vec<GappedHsp>)>,
    posttraceback_lengths: Vec<(usize, usize, i32)>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:72-115;
// c++/src/algo/blast/format/blast_format.cpp:445-477:
// The report of a batch reads these set-up values only; the search returns them
// with or without the comparison-only trace above.
struct LocalBatchSetup {
    parameters: LocalSubjectParameters,
    ungapped_karlin: Vec<KarlinParams>,
    query_validity: Vec<bool>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1011-1024;
// c++/src/algo/blast/format/blast_format.cpp:445-477:
// Formatting consumes the same initial context lengths computed before search.
#[cfg(test)]
pub(super) fn run_local_for_report(
    queries: &[Vec<u8>],
    subjects: &[&[u8]],
    profile: LocalStageDProfile<'_>,
    compo_adjust_mode: BlastCompoAdjustMode,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
) -> Result<(
    Vec<KappaResultHitList>,
    LocalSubjectParameters,
    Vec<KarlinParams>,
    Vec<bool>,
)> {
    run_local_for_report_threads(
        queries,
        SubjectLetters::Records(subjects),
        profile,
        compo_adjust_mode,
        do_sum_stats,
        scoring,
        1,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1414:
// The pool owns one search (or is the run's, `with_run_pool`); ordered subject results enter
// the same Stage D stream.
pub(super) fn run_local_for_report_threads(
    queries: &[Vec<u8>],
    subjects: SubjectLetters<'_>,
    profile: LocalStageDProfile<'_>,
    compo_adjust_mode: BlastCompoAdjustMode,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    threads: usize,
) -> Result<(
    Vec<KappaResultHitList>,
    LocalSubjectParameters,
    Vec<KarlinParams>,
    Vec<bool>,
)> {
    with_search_pool(threads, "tblastn", |pool| {
        let mut setup = None;
        let results = run_local_search_with_pool(
            queries,
            subjects,
            profile,
            compo_adjust_mode,
            do_sum_stats,
            scoring,
            &mut setup,
            None,
            Some(pool),
        )?;
        let setup: LocalBatchSetup = setup.context("missing TBLASTN initial parameters")?;
        Ok((
            results,
            setup.parameters,
            setup.ungapped_karlin,
            setup.query_validity,
        ))
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:277-291
// One query batch searched on a pool that the caller keeps (the port of CLocalBlast::Run for the
// batch, with its report data).
/// `run_local_for_report_threads` on a pool that the caller keeps for all its query batches.
pub(super) fn x_run_local_for_report_in_pool(
    queries: &[Vec<u8>],
    subjects: SubjectLetters<'_>,
    profile: LocalStageDProfile<'_>,
    compo_adjust_mode: BlastCompoAdjustMode,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    pool: &SearchPool<'_>,
) -> Result<XBatchReport> {
    let mut setup = None;
    let results = run_local_search_with_pool(
        queries,
        subjects,
        profile,
        compo_adjust_mode,
        do_sum_stats,
        scoring,
        &mut setup,
        None,
        Some(pool),
    )?;
    let setup: LocalBatchSetup = setup.context("missing TBLASTN initial parameters")?;
    Ok((
        results,
        setup.parameters,
        setup.ungapped_karlin,
        setup.query_validity,
    ))
}

// No NCBI counterpart: names the tuple that run_local_for_report_threads returns for one batch; it
// does not change any value NCBI computes.
/// What `run_local_for_report_threads` returns for one query batch.
pub(super) type XBatchReport = (
    Vec<KappaResultHitList>,
    LocalSubjectParameters,
    Vec<KarlinParams>,
    Vec<bool>,
);

// No NCBI counterpart: reads the LOSAT_X_TBNBATCHRULE variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNBATCHRULE): a switch of the strict set (`utils/xswitch.rs`), on unless
/// it is `0`. It acts only where query batches are
/// searched side by side (`x_side_by_side_batches`).
#[cfg(feature = "parallel")]
fn x_tbn_batch_rule() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNBATCHRULE"))
}

/// EXPERIMENT (LOSAT_X_TBNBATCHRULE): subject bases up to which query batches are searched side
/// by side whatever the number of subjects. Side by side, the serial parts
/// of each batch (query masking, lookup table) overlap; one after another, the pool shares each
/// batch's subjects out. Timed on the workstation (perf-si-20261012, session 17) against
/// the batches one after another (`LOSAT_X_TBNBATCH=0`, a variable that graduated 2026-10-22)
/// with batches a multiple of the threads, at 8 and 16 threads: side by
/// side took 0.26-0.96 of the time on sets of 0.3 to 15.7 Mb (one to 5,000 subjects), 0.94-1.00
/// on 14,000 16S records (20.2 Mb), and 1.00-1.18 on 28,036 16S records (40.6 Mb) and on the
/// fly transcriptome (94 Mb).
#[cfg(feature = "parallel")]
const X_SIDE_BY_SIDE_MAX_BASES: usize = 24_000_000;

// No NCBI counterpart: chooses which query batches are searched side by side, each on one
// thread; every batch is searched as NCBI searches it and the reports keep the batch order.
/// EXPERIMENT (LOSAT_X_TBNBATCHRULE): how many of the first `batches` query batches
/// are searched side by side; the rest are searched one after another with the
/// pool inside each.
///
/// Side by side only when the subjects hold at most `X_SIDE_BY_SIDE_MAX_BASES`, or when there
/// are fewer than two subjects per thread: one after another, each subject is finished (merge,
/// linking, reaping after its last unit) on one thread, so a batch of a few long subjects ends
/// with threads idle. Three fly chromosomes (76.9 Mb) took 1.29 (16 threads) and 1.16 (8) times
/// the time of side by side one after another; the same bases in 11 pieces 1.05 and 0.99, in 21
/// pieces 1.04 and 0.97, in 79 pieces 0.95 and 0.95. Then the batches of the full waves, a
/// multiple of `threads`, since a last wave of fewer batches than threads leaves threads idle
/// until its batches end (10 batches at 8 threads on 16S took 1.53 times as long side by side
/// as one after another; on one 0.6 Mb genome this split took 0.77 of the time of all ten side
/// by side).
#[cfg(feature = "parallel")]
fn x_batches_side_by_side(
    subject_bases: usize,
    subjects: usize,
    batches: usize,
    threads: usize,
) -> usize {
    let many_subjects = subjects >= threads.saturating_mul(2);
    if threads < 2 || (subject_bases > X_SIDE_BY_SIDE_MAX_BASES && many_subjects) {
        return 0;
    }
    batches / threads * threads
}

// No NCBI counterpart: chooses which query batches are searched side by side, each on one
// thread; every batch is searched as NCBI searches it and the reports keep the batch order.
/// Graduated 2026-10-22 (was `LOSAT_X_TBNBATCH`): how many of the first `batches` query
/// batches are searched side by side at `threads` threads: none unless `allowed`
/// (four or more batches, and subject plans that would fit,
/// `XSubjectSide::all_fit`); with `LOSAT_X_TBNBATCHRULE` those of `x_batches_side_by_side`;
/// without it every batch, once there is at least one per thread.
#[cfg(feature = "parallel")]
fn x_side_by_side_batches(
    subject_bases: usize,
    subjects: usize,
    batches: usize,
    threads: usize,
    allowed: bool,
) -> usize {
    if !allowed {
        0
    } else if x_tbn_batch_rule() {
        x_batches_side_by_side(subject_bases, subjects, batches, threads)
    } else if batches >= threads {
        batches
    } else {
        0
    }
}

// No NCBI counterpart: chooses which query batches are searched at the same time on pools of
// their own; every batch is searched as NCBI searches it and the reports keep the batch order.
/// EXPERIMENT (LOSAT_X_TBNOWNPOOLS): how many of the `left` batches that `x_side_by_side_batches`
/// leaves are searched at the same time on pools of their own (`x_batches_on_own_pools`, with
/// `LOSAT_X_TBNOWNPOOLS` and `LOSAT_X_TBNBATCHRULE`; none without them).
#[cfg(feature = "parallel")]
fn x_own_pool_batches(left: usize, subject_bases: usize, subjects: usize, threads: usize) -> usize {
    if x_tbn_own_pools() && x_tbn_batch_rule() {
        x_batches_on_own_pools(left, subject_bases, subjects, threads)
    } else {
        0
    }
}

// No NCBI counterpart: decides whether the subject plans of a search are kept; it does not change
// any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNKEEPFEW): true when a search of `batches` query batches (four or more,
/// as the plans are kept) at `threads` threads searches some of them side by side
/// (`x_side_by_side_batches`) or at the same time on pools of their own (`x_own_pool_batches`).
/// Such a search keeps its plans, as before: batches searched at the same time share one plan of
/// each subject, and those after the first ones read them instead of making them again. Side by
/// side needs plans that fit (at most two bytes per base of `X_SIDE_BY_SIDE_MAX_BASES`, or fewer
/// than two subjects per thread), and own pools fewer than two subjects per thread.
/// Timed on the workstation (2026-10-19, 7 batches of fly proteins against 5,000 fly cDNAs, 14.3
/// Mb, `-db`): at 4 threads (4 batches side by side, then 3 one after another) the search without
/// them took 1.004 of the time (1.003-1.068, 8 pairs); at 8 threads (the 7 batches one after
/// another) 0.998-1.000.
fn x_searched_side_by_side(subjects: SubjectLetters<'_>, batches: usize, threads: usize) -> bool {
    #[cfg(feature = "parallel")]
    {
        if threads < 2 {
            return false;
        }
        let bases = subjects.total_bases();
        let fits = super::search_seed::XSubjectSide::all_fit(bases, subjects.count());
        let side_by_side = x_side_by_side_batches(bases, subjects.count(), batches, threads, fits);
        side_by_side > 0
            || x_own_pool_batches(batches - side_by_side, bases, subjects.count(), threads) > 0
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = (subjects, batches, threads);
        false
    }
}

// No NCBI counterpart: reads the LOSAT_X_TBNOWNPOOLS variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNOWNPOOLS): a switch of the strict set (`utils/xswitch.rs`), on unless
/// it is `0`. It acts only together with
/// `LOSAT_X_TBNBATCHRULE`, on the batches that are not searched side by side, and only natively: in
/// WebAssembly only the main thread can start threads (the WASI thread hosts give a thread a
/// `thread-spawn` that fails, `tests/wasi_thread_host.js` and
/// `web/packages/losat/src/infra/reactor/instance.ts`), and each own pool is built on a thread of its own, so
/// its threads could not start there (os error 6).
#[cfg(feature = "parallel")]
fn x_tbn_own_pools() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    !cfg!(target_family = "wasm")
        && *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNOWNPOOLS"))
}

// No NCBI counterpart: chooses which query batches are searched at the same time on pools of
// their own; every batch is searched as NCBI searches it and the reports keep the batch order.
/// EXPERIMENT (LOSAT_X_TBNOWNPOOLS): how many of the `left` batches that
/// `x_batches_side_by_side` leaves to be searched one after another with the pool inside each
/// are searched instead at the same time, each on a pool of its own with its share of the
/// `threads` (all of them or none).
///
/// Only for fewer than two subjects per thread, and for two to `threads - 1` batches when the
/// subjects hold at most `X_SIDE_BY_SIDE_MAX_BASES` (a pool may then have one thread, as each
/// side-by-side batch has), otherwise for two to `threads / 2` (every pool two threads or more).
/// A batch on a few subjects leaves most of a large pool idle: one 4.7 Mb genome x 16 batches
/// at 32 threads took 0.20 of the time of the batches one after another, x 5 batches at 16
/// threads 0.45 and at 8 threads 0.97, one 0.6 Mb genome x 16 batches at 32 threads 0.15 and x 5
/// batches at 8 threads 0.61, three fly chromosomes (76.9 Mb) x 16 batches at 32 threads 0.81
/// (workstation, perf-si-20261012 session 19). A batch on many subjects keeps a large pool busy,
/// so pools of its own gained nothing there (5,000 16S records x 5 batches at 16 threads 1.09,
/// 28,036 at 16 threads 0.98, the fly transcriptome 1.16), and a pool of one thread is slow
/// where the subjects are long (three fly chromosomes x 5 batches at 8 threads, pools of 2, 2,
/// 2, 1 and 1 threads: 1.44).
#[cfg(feature = "parallel")]
fn x_batches_on_own_pools(
    left: usize,
    subject_bases: usize,
    subjects: usize,
    threads: usize,
) -> usize {
    let few_subjects = subjects < threads.saturating_mul(2);
    let most = if subject_bases <= X_SIDE_BY_SIDE_MAX_BASES {
        threads.saturating_sub(1)
    } else {
        threads / 2
    };
    if few_subjects && left >= 2 && left <= most {
        left
    } else {
        0
    }
}

/// EXPERIMENT (LOSAT_X_TBNOWNPOOLS): the threads of the `index`-th of `pools` pools that share
/// `threads` (the first `threads % pools` get one more).
#[cfg(feature = "parallel")]
fn x_own_pool_threads(threads: usize, pools: usize, index: usize) -> usize {
    threads / pools + usize::from(index < threads % pools)
}

/// Graduated 2026-10-22 (was `LOSAT_X_TBNSSIDE`, `LOSAT_X_TBNBATCH`): the query batches
/// `ranges` of `search` (`tblastn_app.rs`) searched ahead, each with its range, in the order
/// of `ranges`; `None` when `ranges` is empty; the caller's loop then searches its batches
/// (tblastn_app.rs passes none when the database's sequences could not be read).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-289
/// NCBI searches each batch on its own (`CLocalBlast` per batch); a batch reads nothing
/// another batch writes, so only the order of the reports matters, and the caller appends
/// them in the order of `ranges`. One search pool serves every
/// batch instead of one pool per batch; with four or more batches the
/// subject plans (translations, masks) are made once for the search and shared through the
/// pool (a search of fewer batches keeps none, `XSubjectSide::KEEP_FROM_BATCHES`, and with
/// `LOSAT_X_TBNKEEPFEW` nor does one that searches its batches one after another and whose
/// plans exceed what its threads have in flight, `XSubjectSide::kept_within_threads`,
/// `x_searched_side_by_side`); with four or more batches, plans that would fit, and at least
/// one batch per thread the batches are searched side by side, each as the one-thread search. The
/// first failing batch in the order of `ranges` gives the error, as in the reference loop.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_batches<'p>(
    queries: &[Vec<u8>],
    subjects: SubjectLetters<'_>,
    ranges: Vec<std::ops::Range<usize>>,
    profile: impl Fn(std::ops::Range<usize>) -> LocalStageDProfile<'p> + Sync,
    compo_adjust_mode: BlastCompoAdjustMode,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    threads: usize,
) -> Result<Option<Vec<(std::ops::Range<usize>, XBatchReport)>>> {
    let Some(first) = ranges.first() else {
        return Ok(None);
    };
    // The subject plans depend only on the genetic code and the lowercase masking, which
    // every batch's profile shares.
    let first = profile(first.clone());
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:627-631,912-913
    // NCBI frees the translation of a subject when the subject has been searched. The plans are
    // kept only for a later query batch to read, and only from four batches on
    // (XSubjectSide::KEEP_FROM_BATCHES): a search of fewer batches gets no XSubjectSide
    // (XSubjectSide::for_batches) and each plan is dropped after its subject, as in NCBI. With
    // LOSAT_X_TBNKEEPFEW a search of more batches keeps none either unless its plans hold no
    // more than its threads could have in flight (XSubjectSide::kept_within_threads) or it
    // searches batches side by side (x_searched_side_by_side). The number of batches does not
    // decide which loop searches them, and the layout of the batches (side by side, own pools)
    // is decided by the number of batches, whether or not the plans are kept.
    #[cfg(feature = "parallel")]
    let x_ssided_layout = ranges.len() >= super::search_seed::XSubjectSide::KEEP_FROM_BATCHES;
    let x_keep = ranges.len() < super::search_seed::XSubjectSide::KEEP_FROM_BATCHES
        || !super::search_seed::XSubjectSide::keep_few()
        || {
            let (longest, _) = subjects.longest_and_total();
            super::search_seed::XSubjectSide::kept_within_threads(
                longest,
                subjects.total_bases(),
                threads,
            )
        }
        || x_searched_side_by_side(subjects, ranges.len(), threads);
    // The plans are dropped (the switch's own path) or kept as before.
    if ranges.len() >= super::search_seed::XSubjectSide::KEEP_FROM_BATCHES
        && super::search_seed::XSubjectSide::keep_few()
    {
        crate::utils::xreach::x_reach!("LOSAT_X_TBNKEEPFEW", !x_keep);
    }
    let x_subject_side = super::search_seed::XSubjectSide::for_batches(
        x_keep,
        ranges.len(),
        subjects.count(),
        first.genetic_code,
        first.mask_lowercase,
    );
    let search_in = |range: std::ops::Range<usize>, pool: &SearchPool<'_>| {
        x_run_local_for_report_in_pool(
            &queries[range.clone()],
            subjects,
            profile(range),
            compo_adjust_mode,
            do_sum_stats,
            scoring,
            pool,
        )
    };
    // EXPERIMENT (LOSAT_X_TBNOWNPOOLS): the last batches, searched after the shared pool at the
    // same time, each on its own pool (x_batches_on_own_pools).
    #[cfg(feature = "parallel")]
    let mut x_own_pools: Vec<std::ops::Range<usize>> = Vec::new();
    let searched = with_search_pool(threads, "tblastn", |pool| {
        let x_pool_with_plans;
        let pool = match x_subject_side.as_ref() {
            Some(side) => {
                x_pool_with_plans = pool.x_with_ext(side);
                &x_pool_with_plans
            }
            None => pool,
        };
        let mut searched = Vec::with_capacity(ranges.len());
        // Graduated 2026-10-22 (was `LOSAT_X_TBNBATCH`): search the batches side by side, each as
        // the one-thread search. The serial parts of a batch (query masking, lookup table,
        // linking, the hand-over between stages) then no longer leave the other threads
        // idle. It needs four or more batches, and subject plans that
        // would fit (XSubjectSide::all_fit); and at least one batch per thread (with fewer,
        // each batch uses the pool inside). The batches share the kept plans
        // (LOSAT_X_TBNKEEPFEW keeps them for such a search, x_searched_side_by_side); without
        // them each batch would make the plan of the subject it is on and drop it after the
        // subject, as a one-thread search does.
        #[cfg(feature = "parallel")]
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:380-393
        // Dispatch point: NCBI has a mode that searches several query batches at once as nodes of a
        // master node (-mt_mode 1). Here, with four or more batches, plans that
        // would fit and at least one batch per thread, the batches are searched side by side, each
        // as the one-thread search, and the reports are returned in batch order.
        if pool.enabled() {
            let bases = subjects.total_bases();
            let fits = super::search_seed::XSubjectSide::all_fit(bases, subjects.count());
            let side_by_side = x_side_by_side_batches(
                bases,
                subjects.count(),
                ranges.len(),
                pool.threads(),
                x_ssided_layout && fits,
            );
            // The sides are recorded here, where the batches are dispatched (the same helpers also
            // decide whether plans are kept, x_searched_side_by_side).
            crate::utils::xreach::x_reach!("LOSAT_X_TBNBATCH", side_by_side > 0);
            if x_ssided_layout && fits && x_tbn_batch_rule() {
                crate::utils::xreach::x_reach!("LOSAT_X_TBNBATCHRULE", side_by_side > 0);
            }
            // EXPERIMENT (LOSAT_X_TBNOWNPOOLS): the batches left after those, at the same time
            // on pools of their own (x_batches_on_own_pools).
            let own_pools = x_own_pool_batches(
                ranges.len() - side_by_side,
                bases,
                subjects.count(),
                pool.threads(),
            );
            if x_tbn_own_pools() && x_tbn_batch_rule() {
                crate::utils::xreach::x_reach!("LOSAT_X_TBNOWNPOOLS", own_pools > 0);
            }
            if crate::utils::threading::diagnostics_enabled() {
                crate::utils::threading::write_diagnostic_line(format!(
                    "[losat-thread-stage] program=tblastn stage=query_batches batches={} threads={} side_by_side={} one_after_another={} own_pools={}\n",
                    ranges.len(),
                    pool.threads(),
                    side_by_side,
                    ranges.len() - side_by_side - own_pools,
                    own_pools
                ));
            }
            let (ranges, rest) = ranges.split_at(side_by_side);
            let (rest, own) = rest.split_at(rest.len() - own_pools);
            if side_by_side > 0 {
                use rayon::prelude::*;
                // LOSAT_X_TBNKEEPFEW keeps the plans of a search with batches side by side
                // (x_searched_side_by_side).
                debug_assert!(x_subject_side.is_some());
                let x_serial = SearchPool::x_serial();
                let x_serial_with_plans;
                let x_serial = match x_subject_side.as_ref() {
                    Some(side) => {
                        x_serial_with_plans = x_serial.x_with_ext(side);
                        &x_serial_with_plans
                    }
                    None => &x_serial,
                };
                for wave in ranges.chunks(pool.threads() * 4) {
                    let computed: Vec<_> = pool.install(|| {
                        wave.par_iter()
                            .with_max_len(1)
                            .map(|range| search_in(range.clone(), x_serial))
                            .collect()
                    });
                    for (range, report) in wave.iter().zip(computed) {
                        searched.push((range.clone(), report?));
                    }
                }
            }
            // EXPERIMENT (LOSAT_X_TBNBATCHRULE): the batches after the full waves, one
            // after another with the pool inside each.
            for range in rest {
                searched.push((range.clone(), search_in(range.clone(), pool)?));
            }
            x_own_pools = own.to_vec();
            return Ok(Some(searched));
        }
        for range in &ranges {
            searched.push((range.clone(), search_in(range.clone(), pool)?));
        }
        Ok(Some(searched))
    })?;
    #[cfg_attr(
        not(feature = "parallel"),
        expect(
            unused_mut,
            reason = "extended only by the own-pool batches (parallel)"
        )
    )]
    let Some(mut searched) = searched
    else {
        return Ok(None);
    };
    // EXPERIMENT (LOSAT_X_TBNOWNPOOLS): each of the last batches on a pool of its own with its
    // share of the threads, all at the same time, once the search on the shared pool has ended
    // (a pool cannot be started from a worker of another, `with_search_pool`). Each batch is the
    // same search as on the shared pool, whose report does not depend on the number of threads;
    // the reports are taken in batch order, and the first failing batch in that order gives the
    // error.
    #[cfg(feature = "parallel")]
    if !x_own_pools.is_empty() {
        let side = x_subject_side.as_ref();
        let search_in = &search_in;
        let pools = x_own_pools.len();
        let reports: Vec<_> = std::thread::scope(|scope| {
            let handles: Vec<_> = x_own_pools
                .iter()
                .enumerate()
                .map(|(index, range)| {
                    let own = x_own_pool_threads(threads, pools, index);
                    scope.spawn(move || {
                        with_search_pool(own, "tblastn", |own_pool| match side {
                            Some(side) => search_in(range.clone(), &own_pool.x_with_ext(side)),
                            None => search_in(range.clone(), own_pool),
                        })
                    })
                })
                .collect();
            handles.into_iter().map(|handle| handle.join()).collect()
        });
        for (range, report) in x_own_pools.iter().zip(reports) {
            match report {
                Ok(report) => searched.push((range.clone(), report?)),
                Err(panic) => std::panic::resume_unwind(panic),
            }
        }
    }
    Ok(Some(searched))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905;
// c++/src/algo/blast/core/blast_kappa.c:3383-3427,3525-3736,2494-2515;
// c++/src/algo/blast/core/blast_traceback.c:1763-1776:

#[allow(dead_code)] // This internal path is compared before the public CLI is enabled.
pub(super) fn run_local_mode2_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[&[u8]],
    profile: LocalStageDProfile<'_>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        true,
        LocalStageDScoring::default(),
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1481-1501,1644-1707:
// else Blast_TracebackFromHSPList(...); Blast_HSPResultsInsertHSPList(...);
#[allow(dead_code)]
pub(super) fn run_local_mode0_no_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[&[u8]],
    profile: LocalStageDProfile<'_>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        BlastCompoAdjustMode::NoCompositionBasedStats,
        false,
        LocalStageDScoring::default(),
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1501;
// c++/src/algo/blast/core/blast_kappa.c:411-427:
// compositionBasedStats selects redo independently of hitParams->do_sum_stats.
#[allow(dead_code)]
pub(super) fn run_local_mode0_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[&[u8]],
    profile: LocalStageDProfile<'_>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        BlastCompoAdjustMode::NoCompositionBasedStats,
        true,
        LocalStageDScoring::default(),
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1501;
// c++/src/algo/blast/core/blast_kappa.c:411-427:
// compositionBasedStats selects redo independently of hitParams->do_sum_stats.
#[allow(dead_code)]
pub(super) fn run_local_mode2_no_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[&[u8]],
    profile: LocalStageDProfile<'_>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        false,
        LocalStageDScoring::default(),
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:457-463;
// core/blast_kappa.c:2372-2384: selected x-drop bits reach both Stage C
// gapped extension and composition redo.
#[allow(dead_code)]
pub(super) fn run_local_mode2_sum_stats_with_scoring(
    queries: &[Vec<u8>],
    subjects: &[&[u8]],
    profile: LocalStageDProfile<'_>,
    scoring: LocalStageDScoring,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        true,
        scoring,
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:903-936;
// core/blast_traceback.c:198-250,1486-1501:
// the selected protein scoring options enter ordinary traceback and statistics.
#[allow(dead_code)]
pub(super) fn run_local_mode0_no_sum_stats_with_scoring(
    queries: &[Vec<u8>],
    subjects: &[&[u8]],
    profile: LocalStageDProfile<'_>,
    scoring: LocalStageDScoring,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        BlastCompoAdjustMode::NoCompositionBasedStats,
        false,
        scoring,
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1481-1501:
// else Blast_TracebackFromHSPList(...);
fn run_local_search(
    queries: &[Vec<u8>],
    subjects: &[&[u8]],
    profile: LocalStageDProfile<'_>,
    compo_adjust_mode: BlastCompoAdjustMode,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    observer: Option<&mut StageDBoundaryTrace>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search_with_pool(
        queries,
        SubjectLetters::Records(subjects),
        profile,
        compo_adjust_mode,
        do_sum_stats,
        scoring,
        &mut None,
        observer,
        None,
    )
}

/// NCBI's engine error of a query batch's setup (`BLAST_GapAlignSetUp`, c++/src/algo/blast/core/blast_setup.c:888-1000):
/// NCBI 2.17.0 tblastn writes the start of the outfmt 0 report (and the reports of the batches
/// before), then `BLAST engine error: <text>`, exit status 3. `tblastn_app::batch_failure` fails
/// the batch with it, as a read of the database that fails. The text is one of the constants.
#[derive(Debug)]
pub(crate) struct SetupEngineError(pub(crate) &'static str);

impl SetupEngineError {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:433-437
    /// `BlastExtensionParametersNew` (blast_setup.c:960, before the subject length) finds no valid
    /// `sbp->kbp`: a PSSM query whose file's ungapped values make `kbp_psi` invalid
    /// (`query_set_setup`; probes `x.ungneg.*`, `r3.*ungneg.*` of the PS7 review).
    pub(crate) const NO_VALID_KARLIN_BLOCK: Self =
        Self(crate::report::query_warnings::INVALID_QUERY_MESSAGE);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-217
    /// `BLASTERR_SUBJECT_LENGTH_INVALID` (`run_local_search_with_pool`, blast_setup.c:969-980):
    /// `-subject` of three empty records, or of `A`, `` and ``, with `-searchsp` or `OLD_FSC`
    /// (env parity probe f_probe of 2026-11-12).
    pub(crate) const SUBJECT_LENGTH_INVALID: Self =
        Self(crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT);

    /// `BlastHitSavingParametersNew`'s status 1 (`stats::setup::hit_saving_parameters_fail`):
    /// sum statistics with a setup subject length that is not positive, the average of a negative
    /// `NSEQ` without the Gumbel block (NCBI 2.17.0 `-searchsp` over `NSEQ 4294967295`, probe
    /// run1 of line RG D1a-3).
    pub(crate) const HIT_SAVING_PARAMETERS: Self =
        Self(crate::stats::setup::HIT_SAVING_PARAMETERS_ERROR);
}

impl std::fmt::Display for SetupEngineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{ENGINE_ERROR} {}", self.0)
    }
}

impl std::error::Error for SetupEngineError {}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-301
/// LOSAT's refusal of a batch's setup where NCBI's crashes (`common::database::check_setup_count`:
/// the statistics' count of 0 without the Gumbel block). NCBI writes the reports of each batch
/// before the next is searched, so its tabular output before the crash holds those of the batches
/// before (NCBI 2.17.0, D1a-5 of line RG: `OLD_FSC=1 BATCH_SIZE=1` over `NSEQ 0`, two queries
/// without letters to search then a valid one, the 152 bytes of their outfmt 7 reports without the
/// last line; the crash cuts the buffered output of the other formats short):
/// `tblastn_app::batch_failure` makes the refusal after them, as an engine error of the batch.
#[derive(Debug)]
pub(crate) struct SetupRefusal(pub(crate) anyhow::Error);

impl std::fmt::Display for SetupRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for SetupRefusal {}

/// Diagnostic (LOSAT_TIMING): adds the time from `start` to its drop to the common line's `search`
/// (utils/xtiming.rs). No NCBI counterpart; it reads the clock only with LOSAT_TIMING.
struct XSearchTimer(Option<std::time::Instant>);

impl XSearchTimer {
    fn start() -> Self {
        Self(crate::utils::xtiming::clock())
    }
}

impl Drop for XSearchTimer {
    fn drop(&mut self) {
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Search, self.0);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1414,1469-1475:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905:
// Preliminary work uses one immutable subject; linking and stream writes follow OID order.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039 (cited at
// `match_redo`)
// Dispatch point: the redo of a match is chosen once per query batch, by NCBI's `eTbackExt`
// (`-use_sw_tback`): the Smith-Waterman redo (`match_redo::SmithWatermanRedo`) or the X-drop redo
// of composition-based statistics (`match_redo::XDropRedo`); each is its own instance of the
// batch's search (`run_local_search_with_redo`).
#[allow(clippy::too_many_arguments)]
fn run_local_search_with_pool(
    queries: &[Vec<u8>],
    subjects: SubjectLetters<'_>,
    profile: LocalStageDProfile<'_>,
    compo_adjust_mode: BlastCompoAdjustMode,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    setup: &mut Option<LocalBatchSetup>,
    observer: Option<&mut StageDBoundaryTrace>,
    pool: Option<&SearchPool<'_>>,
) -> Result<Vec<KappaResultHitList>> {
    if scoring.smith_waterman {
        return run_local_search_with_redo(
            queries,
            subjects,
            profile,
            compo_adjust_mode,
            do_sum_stats,
            scoring,
            setup,
            observer,
            pool,
            &SmithWatermanRedo::new(scoring.matrix, compo_adjust_mode, profile.expect_value),
        );
    }
    run_local_search_with_redo(
        queries,
        subjects,
        profile,
        compo_adjust_mode,
        do_sum_stats,
        scoring,
        setup,
        observer,
        pool,
        &XDropRedo,
    )
}

/// The search of one query batch (`run_local_search_with_pool`) with the redo of a match `R`.
#[allow(clippy::too_many_arguments)]
fn run_local_search_with_redo<R: MatchRedo>(
    queries: &[Vec<u8>],
    subjects: SubjectLetters<'_>,
    profile: LocalStageDProfile<'_>,
    compo_adjust_mode: BlastCompoAdjustMode,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    setup: &mut Option<LocalBatchSetup>,
    mut observer: Option<&mut StageDBoundaryTrace>,
    pool: Option<&SearchPool<'_>>,
    match_redo: &R,
) -> Result<Vec<KappaResultHitList>> {
    ensure!(
        !queries.is_empty(),
        "TBLASTN Stage D requires a protein query"
    );
    // A database search has no subjects in hand when its sequences could not be read: its
    // batches are set up and reported, and `tblastn_app::search` fails the first with a valid
    // context.
    ensure!(
        subjects.count() > 0 || profile.subject_set.database.is_some(),
        "TBLASTN Stage D requires a local subject"
    );
    ensure!(
        profile.max_target_seqs > 0,
        "TBLASTN hitlist size must be positive"
    );
    let code = GeneticCode::try_from_id(profile.genetic_code).map_err(anyhow::Error::msg)?;
    // The set's numbers: a database's, or those of the FASTA subjects (`SubjectSetNumbers`).
    let total_nt = profile.subject_set.total_nt(&subjects);
    let stats_nt = profile.subject_set.eff_length_nt(total_nt);
    let num_subjects = profile.subject_set.num_seqs(&subjects);
    let query_lengths: Vec<i32> = queries
        .iter()
        .map(|query| i32::try_from(query.len()))
        .collect::<std::result::Result<_, _>>()?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:908-936;
    // core/blast_stat.c:2547-2605: reject scoring schemes absent from
    // the pinned NCBI gapped table before any search.
    let spec = ProteinScoringSpec {
        matrix: scoring.matrix,
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
    };
    ensure!(
        protein_scoring_supported(&spec),
        "unsupported TBLASTN protein scoring scheme"
    );
    let gapped = lookup_protein_params(&spec);
    let ungapped = lookup_protein_params_ungapped(scoring.matrix);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4091-4171;
    // core/blast_hits.c:1870-1895: gbp is optional; absent gbp selects
    // BLAST_KarlinStoE_simple with the context's effective search space.
    let gumbel = profile.subject_set.gumbel_block(&spec, total_nt);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_db_adapter.cpp:133-135
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:226-240
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-973
    // The local `-subject` sequences are a multi-sequence source: the shortest record (an
    // empty one counts as 0, in nucleotides) is raised to 10 before the division. A database
    // answers `GetMinLength` (`SubjectSetNumbers::shortest_subject_length`); a quotient of 0 is
    // NCBI's engine error of the first searched batch (`tblastn_app.rs`).
    // Without a Gumbel block it is the average, in nucleotides (blast_setup.c:975). The block
    // itself decides, as `sbp->gbp` does in NCBI, here and in the average-length error below.
    // `None`: the setup fails (`stats::setup::setup_subject_length`, below).
    let min_subject_length = if gumbel.is_none() {
        crate::stats::setup::setup_subject_length(total_nt, num_subjects, None)
    } else {
        crate::stats::setup::setup_subject_length(
            total_nt,
            num_subjects,
            Some((profile.subject_set.shortest_subject_length(&subjects), true)),
        )
    };
    let options_searchsp = profile.subject_set.searchsp;
    let QuerySetSetup {
        working_frames,
        contexts,
        valid_contexts,
        gapped_by_context,
        parameters,
        recipe,
    } = query_set_setup(
        queries,
        QuerySetInputs {
            profile,
            scoring,
            gapped,
            ungapped: &ungapped,
            gumbel: gumbel.as_ref(),
            total_nt: stats_nt,
            num_subjects,
            min_subject_length: min_subject_length.unwrap_or(0),
            compo_adjust_mode,
            do_sum_stats,
        },
        options_searchsp.as_ref().map(std::slice::from_ref),
        true,
        pool,
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
    // Without the Gumbel block (`-searchsp` or `OLD_FSC`; LOSAT's tblastn has no ungapped
    // search, which has none either, blast_setup.c:517-527) a set of subjects whose average is
    // below one letter, in nucleotides, fails the setup of a batch that NCBI searches (one with
    // a valid context): `tblastn_app::search` fails it with NCBI's engine error after the
    // reports of the batches before (`SetupEngineError::SUBJECT_LENGTH_INVALID`). A count of 0
    // (of a positive total) is the division by 0 where NCBI crashes, LOSAT's refusal there
    // (`common::database::check_setup_count`, D1a-5 of line RG), after the same reports
    // (`SetupRefusal`); a run whose batches have no valid context does not reach it and is
    // reported as NCBI's. `query_set_setup` above (NCBI's `BLAST_CalcEffLengths`,
    // blast_setup.c:946-951, before the division at :975) refuses first where NCBI's search
    // space overflows: with a count of 0 as well, that LOSAT refusal takes the place of this one
    // for a run that NCBI does not finish. It is not a batch failure, so it comes before the
    // reports of the batches before, which NCBI writes before it crashes (line RG D1a-6,
    // review A F2).
    if valid_contexts.iter().any(|&valid| valid) {
        crate::algorithm::common::database::check_setup_count(
            total_nt,
            num_subjects,
            gumbel.is_some(),
            "TBLASTN",
        )
        .map_err(SetupRefusal)?;
        match min_subject_length {
            None => return Err(SetupEngineError::SUBJECT_LENGTH_INVALID.into()),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:759-760
            // (cited at `stats::setup::hit_saving_parameters_fail`); TBLASTN is gapped.
            Some(length)
                if crate::stats::setup::hit_saving_parameters_fail(do_sum_stats, true, length) =>
            {
                return Err(SetupEngineError::HIT_SAVING_PARAMETERS.into())
            }
            Some(_) => {}
        }
    }
    let batch_setup = LocalBatchSetup {
        parameters: parameters.clone(),
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:72-115:
        ungapped_karlin: contexts
            .iter()
            .map(|context| context.karlin_params)
            .collect(),
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-97: no valid context returns before
        // Karlin block initialization and leaves search space at zero.
        query_validity: contexts.iter().map(|context| context.is_valid).collect(),
    };
    if let Some(ref mut trace) = observer {
        trace.parameters = Some(batch_setup.parameters.clone());
        trace.ungapped_karlin = batch_setup.ungapped_karlin.clone();
        trace.query_validity = batch_setup.query_validity.clone();
        trace.subject_nt_lengths = subjects.x_read_lengths()?;
    }
    *setup = Some(batch_setup);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180,197-224:
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:322-327:
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2815-2824:
    // Preserve the parameter/validity state for reporting, then skip subject
    // search when every query context has failed Karlin setup.
    if !valid_contexts.iter().any(|&valid| valid) || subjects.count() == 0 {
        return (0..queries.len())
            .map(|_| KappaResultHitList::new(profile.max_target_seqs))
            .collect();
    }
    // Diagnostic (LOSAT_TIMING): from here to the end of the batch's search (the preliminary
    // search of all subjects, which translates them and builds the lookup table inside, the
    // composition redo and the traceback) is `search` of the common line (utils/xtiming.rs); a
    // guard, as the function has many returns. No NCBI counterpart.
    let _x_search_timer = XSearchTimer::start();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-220 (cited at
    // `stats::setup::preliminary_option_space`): the preliminary stage's option spaces, which its
    // setup reads in place of its own where the statistics' total before the fallback to OID 0 is
    // not positive, and which the statistics of each subject read
    // (`SubjectSetNumbers::per_subject`, blast_engine.c:1434-1443).
    let per_subject = profile.subject_set.per_subject();
    let option_spaces: Option<Vec<i64>> = match profile.subject_set.database {
        Some(database) if database.gumbel_db_length <= 0 || per_subject => Some(
            recipe
                .query_contexts
                .iter()
                .zip(&recipe.lengths)
                .map(|(&(query_length, valid), ka)| {
                    crate::stats::setup::preliminary_option_space(
                        if valid { query_length as i64 } else { 0 },
                        database.gumbel_db_length,
                        database.request_count,
                        profile.subject_set.dbsize.unwrap_or(0),
                        options_searchsp.unwrap_or(0),
                        true,
                        ka,
                        "TBLASTN",
                    )
                })
                .collect::<Result<_, _>>()?,
        ),
        _ => None,
    };
    let preliminary_parameters = match (&option_spaces, profile.subject_set.database) {
        (Some(spaces), Some(database)) if database.gumbel_db_length <= 0 => recipe.parameters(
            stats_nt,
            Some(spaces),
            LocalParameterCall::Initial {
                min_subject_length: min_subject_length.unwrap_or(0),
                composition_based_stats: compo_adjust_mode.ncbi_value(),
            },
        )?,
        _ => parameters.clone(),
    };
    let link = preliminary_parameters.link;
    let cutoff_arrays = ContextCutoffArrays::of(&preliminary_parameters);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-458: the gapped X-drop from `kbp_gap` (a PSSM query's
    // `kbp_gap_psi`, `query_set_setup`).
    let search = preliminary_profile(
        profile,
        scoring,
        profile
            .pssm
            .map_or(gapped, |psi| psi.gapped_karlin(&gapped)),
        &cutoff_arrays,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:43-70;
    // c++/src/algo/blast/core/hspfilter_collector.c:105-161:
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:43-70:
    // Computed without NCBI's `int` overflow (`blastn/hsp.rs` `get_prelim_hitlist_size`).
    let prelim_size = crate::algorithm::blastn::hsp::get_prelim_hitlist_size(
        profile.max_target_seqs,
        compo_adjust_mode.is_enabled(),
        true,
    );
    let query_refs: Vec<_> = queries.iter().map(Vec::as_slice).collect();
    let protein_split =
        crate::algorithm::blastp::protein_query_split::protein_query_split_sizes(20000)?;
    let batch_lengths: Vec<usize> = queries.iter().map(Vec::len).collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:122-129 (cited at
    // `blastp::blast_engine::psi_round::run_queries`)
    // A PSSM query (`-in_pssm`: `CPsiBlast` searches the PSSM's query through a
    // `CObjMgrFree_QueryFactory`, psiblast_impl.cpp:115-119) is never split; the sizes above are
    // still read, as NCBI's query splitter reads them first.
    let split = match profile.pssm {
        Some(_) => None,
        None => split_protein_batch(&batch_lengths, protein_split),
    };
    let preliminary_hitlists = match split {
        None => run_preliminary_subjects(
            &PreliminaryInputs {
                queries: &query_refs,
                frames: &working_frames,
                genetic_code: profile.genetic_code,
                search,
                link,
                query_lengths: &query_lengths,
                lengths: &preliminary_parameters.lengths,
                gapped_by_context: &gapped_by_context,
                gumbel: gumbel.as_ref(),
                prelim_evalue: preliminary_parameters.prelim_evalue,
                prelim_size,
                culling_limit: profile.culling_limit,
                pssm: profile.pssm.map(|psi| &psi.pssm),
                update: per_subject.then(|| SubjectUpdate {
                    recipe: &recipe,
                    option_spaces: option_spaces.as_deref().unwrap_or(&[]),
                    dbsize_nt: profile.subject_set.dbsize_nt(),
                    min_length: scoring.word_size * 3 + 2,
                }),
            },
            subjects,
            pool,
            observer.as_deref_mut(),
        )?,
        Some(chunks) => split_preliminary_hitlists(
            queries,
            &chunks,
            protein_split,
            &parameters,
            option_spaces.as_deref(),
            QuerySetInputs {
                profile,
                scoring,
                gapped,
                ungapped: &ungapped,
                gumbel: gumbel.as_ref(),
                total_nt: stats_nt,
                num_subjects,
                min_subject_length: min_subject_length.unwrap_or(0),
                compo_adjust_mode,
                do_sum_stats,
            },
            subjects,
            prelim_size,
            pool,
        )?,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1839-1846
    // `ADAPTIVE_CBS`: the stream's close trims each query's preliminary hit list
    // (`core::blast_hspstream`).
    let mut preliminary_hitlists = preliminary_hitlists;
    if crate::core::blast_hspstream::adaptive_cbs_close(
        compo_adjust_mode.is_enabled(),
        profile.max_target_seqs,
    ) {
        for hitlist in preliminary_hitlists.iter_mut().flatten() {
            hitlist.cbs_close_trim(profile.max_target_seqs);
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:144-152,289-319:
    // Blast_HSPResultsReverseSort(results); stream reads from the end.
    let mut preliminary_lists: Vec<_> = preliminary_hitlists
        .iter()
        .enumerate()
        .flat_map(|(context, hitlist)| {
            hitlist.iter().flat_map(move |hitlist| {
                hitlist
                    .lists()
                    .iter()
                    .map(move |list| (context, list.oid as usize, list.hsps.clone()))
            })
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3346-3357
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:657-660
    // The redo's set-up reads the frequency ratios of the matrix (`s_GetAlignParams`, its
    // `s_MatrixInfoInit`), which IDENTITY has none of (`_PSIMatrixFrequencyRatiosNew` returns NULL
    // past its named matrices, matrix_freq_ratios.c:1734-1736): the batch's traceback fails before
    // the redo reads a match, whatever the batch found (`Traceback failed`; no heap holds a match).
    // Only `-use_sw_tback` gets here with IDENTITY, which NCBI searches in mode 0
    // (`identity_cbs_reset`), where the X-drop path runs no redo.
    if (compo_adjust_mode.is_enabled() || R::SMITH_WATERMAN)
        && compo_matrix_data(scoring.matrix).is_none()
    {
        return Err(sw_batch_failed(false));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96,158-204,
    // 569-610: the ordinary traceback stream sorts HSP lists by decreasing
    // OID, then each batch reads from the end; equal-OID query lists were
    // inserted in query order and are consumed in reverse query order.
    // Composition redo uses the score-sorted stream at lines 289-319.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:674-681 (the score sort
    // only with composition-based statistics: the Smith-Waterman redo of `-comp_based_stats 0`
    // reads the stream in this OID order too)
    if compo_adjust_mode.is_enabled() {
        preliminary_lists.sort_by(|(qa, oa, a), (qb, ob, b)| {
            qa.cmp(qb)
                .then_with(|| compare_preliminary_lists_for_kappa(*oa as i32, a, *ob as i32, b))
        });
    } else {
        preliminary_lists.sort_by(|(qa, oa, _), (qb, ob, _)| oa.cmp(ob).then_with(|| qb.cmp(qa)));
    }

    // The subject of a list in the stream, for NCBI's fetches with `check_oid_exclusion`.
    let list_ordinal = |(_, oid, _): &(usize, usize, LinkedHspList)| u32::try_from(*oid).ok();
    let fetch_refusal = |failed| FetchExcluded::refusal(failed, "TBLASTN");
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1498 (cited at
    // `match_redo`): the Smith-Waterman redo runs without composition-based statistics too.
    if !compo_adjust_mode.is_enabled() && !R::SMITH_WATERMAN {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1567-1577 (cited at
        // `blastp::blast_engine::database::exclude_fetched`): the traceback's fetches (`check_oid_exclusion`) exclude a subject that
        // the negative seqid list keeps without a line; its lists leave here, after the trim of
        // `ADAPTIVE_CBS` (tblastn is gapped: `-ungapped` is refused).
        if let Some(excluded) = profile.fetch_excluded {
            excluded
                .retain(&mut preliminary_lists, list_ordinal)
                .map_err(fetch_refusal)?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1481-1501
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1644-1707
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1776
        // The control path performs ordinary traceback and posttraceback statistics, then result
        // post-pipes.
        return finish_local_mode0_no_sum_stats(
            queries,
            subjects,
            profile,
            preliminary_lists,
            &query_lengths,
            &gapped_by_context,
            gumbel.as_ref(),
            &code,
            &parameters.lengths,
            link.as_ref(),
            scoring,
            observer.as_deref_mut(),
            pool,
        );
    }
    // The redo's fetches (`s_MatchingSequenceInitialize`, `check_oid_exclusion`) exclude the
    // same subjects; their lists stay in the stream as markers (`FetchMarks`): the loops below test
    // a marker's early termination as any list's, and one that passes clears NCBI's status of the
    // Smith-Waterman redo (`SwBatchStatus::excluded_at_fetch`) and is neither redone nor kept.
    let excluded_at_fetch =
        FetchMarks::of(profile.fetch_excluded, &preliminary_lists, list_ordinal)
            .map_err(fetch_refusal)?;
    // A PSSM query's Smith-Waterman redo (`-in_pssm`, composition-based statistics 0) is refused at
    // its first match, the first list that the fetch keeps (the heaps are empty there, so early
    // termination lets it go on): `match_redo::sw_pssm_refusal`.
    if R::SMITH_WATERMAN
        && profile.pssm.is_some()
        && (0..preliminary_lists.len()).any(|index| !excluded_at_fetch.excluded(index))
    {
        return Err(sw_pssm_refusal());
    }
    // With a `-db` stream, the subjects of the lists that the redo fetches, each read
    // when first aligned, held until their last match is done (LOSAT_X_TBNSTAGEDROP; with it off
    // until the redo ends).
    let x_stage = super::x_db_stream::XStageSubjects::new(
        subjects,
        preliminary_lists
            .iter()
            .enumerate()
            .filter(|(index, _)| !excluded_at_fetch.excluded(*index))
            .map(|(_, (_, oid, _))| *oid),
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2352-2384,2418-2479,
    // 3099-3110,3290-3302:
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-624;
    // c++/src/algo/blast/core/blast_kappa.c:2309-2337,3244-3248,515-526:
    // For TBLASTN s_ComputeNumIdentities passes the working query sequence.
    // Stage C uses these same hard-masked bytes for extension and traceback.
    let encoded: Vec<Vec<u8>> = working_frames
        .iter()
        .map(|frames| frames[0].aa_seq[1..1 + frames[0].aa_len].to_vec())
        .collect();
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:983-1021:
    // The formatter reads the original query, not the hard-masked Kappa buffer.
    let report_queries: Vec<Vec<u8>> = queries
        .iter()
        .map(|query| {
            query
                .iter()
                .copied()
                .map(|aa| aa_char_to_ncbistdaa(aa.to_ascii_uppercase()))
                .collect()
        })
        .collect();
    let query_infos: Vec<_> = encoded
        .iter()
        .enumerate()
        .map(|(index, seq)| BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(seq),
            composition: read_aa_composition(seq),
            eff_search_space: parameters.lengths[index].eff_searchsp as f64,
            words: Some(build_query_word_hashes(seq)),
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2428,2450-2479:
    // do_link_hsps and cutoff_s follow the sum statistics of the hit-saving parameters
    // (`hitParams->do_sum_stats`: the search's linking parameters, `parameters.link`).
    let gumbel = gumbel.as_ref();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039 (cited at
    // `match_redo`): the redo's `smithWaterman`.
    let redo_params = BlastRedoAlignParams {
        smith_waterman: R::SMITH_WATERMAN,
        ..local_kappa_redo_params(
            scoring.matrix,
            scoring.gap_open,
            scoring.gap_extend,
            &gapped_by_context,
            &valid_contexts,
            &parameters,
            *query_lengths.iter().max().unwrap(),
            compo_adjust_mode,
            false,
            scoring.final_xdrop_bits,
            local_extension_final_xdrop(
                scoring.gap_xdrop_bits,
                scoring.final_xdrop_bits,
                gapped.lambda,
            )?,
        )?
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3108
    // `localScalingFactor`, the redo parameters' (32 with
    // composition-based statistics, 1 without them, where only the Smith-Waterman redo runs).
    let scale = redo_params.score_divisor;
    let scaled = KarlinParams {
        lambda: gapped.lambda / scale,
        ..gapped
    };
    // `kbp->logK` of the redo's block (the scaling leaves K): read by the Smith-Waterman redo of
    // e-value mode and by the normalization of each kept match.
    let log_k = gapped.k.glibc_ln();
    let scaled_by_context = vec![scaled; queries.len()];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:411-427;
    // core/blast_parameters.c:1245-1253: postredo link, when enabled, uses
    // the updated zero small-gap cutoff; direct E-value uses raw subject length.
    let mut post_link = link;
    if let Some(link) = post_link.as_mut() {
        link.cutoff_small_gap = 0;
    }
    let mut heaps: Vec<_> = (0..queries.len())
        .map(|_| CompoHeap::new(profile.max_target_seqs, 0.002))
        .collect::<Result<_>>()?;
    let mut owned_lists = HashMap::new();
    let mut scratch = match_redo.scratch();
    let mut workspace = BlastCompositionWorkspace::new_blosum62();
    if scoring.matrix != ScoringMatrix::Blosum62 {
        use_matrix_workspace(&mut workspace, scoring.matrix)?;
    }
    // Everything one (query, subject) match needs
    // between the early-termination test and the heap insertion, as a function of
    // the match alone.  NCBI c++/src/algo/blast/core/blast_kappa.c:3525-3736 keeps
    // the same state thread-local (matrix, NRrecord, gapAlign) for its redo loop.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3534,3625-3626
    // redo_item is the redo of one (query, subject) match between the early-termination test and
    // the heap insertion, taken out as a function of the match. The steps are those of the serial
    // loop below.
    let redo_item = |context: usize,
                     oid: usize,
                     list: &LinkedHspList,
                     redo_params: &BlastRedoAlignParams,
                     scratch: &mut R::Scratch,
                     workspace: &mut BlastCompositionWorkspace,
                     heaps: Option<&mut [CompoHeap]>|
     -> Result<Option<(CompoHeapRecord, KappaResultList)>> {
        // The subject's letters, or (EXPERIMENT, LOSAT_X_TBNSTORE) its packed form: the redo and
        // the translations below read windows of it, as NCBI's do. A
        // `-db` sequence is read once for the redo (`x_stage`); EXPERIMENT (LOSAT_X_TBNSTAGEDROP):
        // the stage lets it go after its last match (`_done`, declared before `held`).
        let _done = x_stage.x_match(oid);
        let held = x_stage.held(oid)?;
        let subject = held.nt();
        let incoming: Vec<_> = list.hsps.iter().map(|hsp| hsp.hsp).collect();
        // The matrix that the match's adjustments leave, for the engine positives
        // (`postredo_engine_positives`).
        let mut matrix_state = None;
        let mut redo = match_redo.redo(
            &incoming,
            i32::try_from(context)?,
            &query_infos,
            subject,
            profile.genetic_code,
            redo_params,
            scaled.lambda,
            log_k,
            scoring.matrix,
            oid,
            heaps,
            scratch,
            workspace,
            &mut matrix_state,
        )?;
        let converted = convert_distinct_alignments(&mut redo.alignments_by_query[context])?;
        if converted.is_empty() {
            return Ok(None);
        }
        let mut hsp_values: Vec<_> = converted.iter().map(|h| (h.context, h.hsp)).collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3661-3706:
        let survivor_indices = reap_contained_postredo_hsps(&mut hsp_values);
        let mut slots: Vec<_> = converted.into_iter().map(Some).collect();
        let survivors: Vec<_> = survivor_indices
            .into_iter()
            .map(|index| slots[index].take().expect("unique containment survivor"))
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:411-427,364-368:
        // else Blast_HSPListGetEvalues(..., s_GetSubjectLength(...), ..., 1.0);
        let mut postredo = if let Some(link) = post_link.as_ref() {
            link_preliminary_hsps(
                &hsp_values,
                &query_lengths,
                &parameters.lengths,
                i32::try_from(subject.len())?,
                &scaled_by_context,
                gumbel,
                link,
            )?
        } else {
            direct_local_evalues(
                &hsp_values,
                &query_lengths,
                i32::try_from(subject.len())?,
                &scaled_by_context,
                &parameters.lengths,
                gumbel,
            )
        };
        reap_by_evalue(&mut postredo, profile.expect_value);
        if postredo.hsps.is_empty() {
            return Ok(None);
        }
        let candidate = CompoHeapRecord {
            subject_index: i32::try_from(oid)?,
            best_evalue: postredo.best_evalue,
            best_score: postredo.hsps[0].hsp.score,
        };
        let bits = normalize_postredo_scores(&mut postredo, scaled.lambda, log_k, scale);
        let mut target = TargetTranslation::of(subject, &code);
        let mut owned: Vec<_> = survivors.into_iter().map(Some).collect();
        let payloads = postredo
            .hsps
            .iter()
            .enumerate()
            .map(|(index, hsp)| {
                let converted = owned[hsp.source_index]
                    .take()
                    .context("TBLASTN Kappa source owner reused")?;
                let (ident, positives, length, mismatches, gap_opens, gap_letters) =
                    postredo_converted_stats_with_matrix(
                        &converted,
                        &encoded[context],
                        &mut target,
                        scoring.matrix,
                    )?;
                // With the Seq-aligns requested, NCBI's engine positives: by the matrix of the
                // match's last adjustment (`postredo_engine_positives`).
                let positives = match matrix_state.as_ref().filter(|_| profile.engine_positives) {
                    Some(adjusted) => postredo_engine_positives(
                        &converted,
                        &encoded[context],
                        &mut target,
                        adjusted,
                    )?,
                    None => positives,
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:515-526 uses the masked working query;
                // objtools/align_format/tabular.cpp:983-1021 recomputes the
                // report counts from the original aligned query sequence.
                let (report_ident, report_positives, _, report_mismatches, _, _) =
                    postredo_converted_stats_with_matrix(
                        &converted,
                        &report_queries[context],
                        &mut target,
                        scoring.matrix,
                    )?;
                Ok(KappaHspPayload {
                    bit_score: bits[index],
                    num_ident: i32::try_from(ident)?,
                    num_positives: positives,
                    report_num_ident: report_ident,
                    report_num_positives: report_positives,
                    report_mismatches: i64::try_from(report_mismatches)?,
                    align_length: length,
                    mismatches,
                    gap_opens,
                    gap_letters,
                    edit_script: converted.edit_script,
                    matrix_adjust_rule: converted.matrix_adjust_rule,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Some((
            candidate,
            KappaResultList {
                oid: candidate.subject_index,
                hsps: postredo,
                payloads,
            },
        )))
    };
    // The heap update of one redone match, in the order of the preliminary lists.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3705-3715
    // The heap update of one redone match (BlastCompo_HeapWouldInsert / HeapInsert), applied in the
    // order of the preliminary lists.
    let mut keep_item = |context: usize,
                         heaps: &mut Vec<CompoHeap>,
                         item: Option<(CompoHeapRecord, KappaResultList)>| {
        let Some((candidate, result)) = item else {
            return;
        };
        if heaps[context].would_insert(candidate) {
            let discarded = heaps[context].insert(candidate);
            owned_lists.insert((context, candidate.subject_index), result);
            if let Some(previous) = discarded {
                owned_lists.remove(&(context, previous.subject_index));
            }
        }
    };
    // Redo the matches of one batch on the search
    // pool, then apply early termination and the heap updates in list order.  A
    // match is a function of its own preliminary list; the only state shared
    // between matches is the heap array, which only this thread touches.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3448-3451
    // NCBI redoes the matches in this omp loop. The matches of a batch are
    // redone on the pool by redo_item and then applied in list order (early termination is tested
    // again in order). A match depends only on its own preliminary list.
    // NCBI's status of the last match that the Smith-Waterman redo reached (`status_code_tld`,
    // `match_redo::match_outcome`); the X-drop redo has none.
    let mut sw_status = SwBatchStatus::default();
    // The Smith-Waterman redo of e-value mode reads the heaps (`match_redo::MatchRedo::
    // reads_heaps`): its matches are redone in stream order on this thread, as NCBI's
    // `-num_threads 1` redoes them.
    #[cfg(feature = "parallel")]
    let parallel_redo = pool.is_some_and(SearchPool::enabled)
        && preliminary_lists.len() > 1
        && !match_redo.reads_heaps(&redo_params);
    #[cfg(not(feature = "parallel"))]
    let parallel_redo = false;
    // No NCBI counterpart: stage report for the thread-use statistics; it does not change any value
    // NCBI computes.
    crate::utils::threading::report_stage(
        "tblastn",
        "kappa_redo",
        preliminary_lists.len(),
        parallel_redo,
    );
    if parallel_redo {
        #[cfg(feature = "parallel")]
        {
            let pool = pool.expect("parallel pool selected");
            let make_params = || {
                local_kappa_redo_params(
                    scoring.matrix,
                    scoring.gap_open,
                    scoring.gap_extend,
                    &gapped_by_context,
                    &valid_contexts,
                    &parameters,
                    *query_lengths.iter().max().unwrap(),
                    compo_adjust_mode,
                    false,
                    scoring.final_xdrop_bits,
                    local_extension_final_xdrop(
                        scoring.gap_xdrop_bits,
                        scoring.final_xdrop_bits,
                        gapped.lambda,
                    )?,
                )
            };
            let states = (0..pool.threads())
                .map(|_| {
                    let mut worker = XRedoWorker {
                        params: BlastRedoAlignParams {
                            smith_waterman: R::SMITH_WATERMAN,
                            ..make_params()?
                        },
                        scratch: match_redo.scratch(),
                        workspace: BlastCompositionWorkspace::new_blosum62(),
                    };
                    if scoring.matrix != ScoringMatrix::Blosum62 {
                        use_matrix_workspace(&mut worker.workspace, scoring.matrix)?;
                    }
                    Ok(std::sync::Mutex::new(worker))
                })
                .collect::<Result<Vec<_>>>()?;
            let batch_size = pool.threads().saturating_mul(32).max(64);
            for (chunk, batch) in preliminary_lists.chunks(batch_size).enumerate() {
                let first = chunk * batch_size;
                // Early termination only becomes true as the heaps fill, so a match
                // that is skipped now is also skipped at its turn; anything else is
                // rechecked in order below. A marker (`excluded_at_fetch`) is never redone.
                let wanted: Vec<bool> = batch
                    .iter()
                    .enumerate()
                    .map(|(at, (_, _, list))| {
                        !excluded_at_fetch.excluded(first + at)
                            && !compo_early_termination(list.best_evalue, &heaps)
                    })
                    .collect();
                let lpt = x_tbn_redo_lpt() && {
                    let acts = wanted.iter().filter(|&&wanted| wanted).count() > 1;
                    crate::utils::xreach::x_reach!("LOSAT_X_TBNREDOLPT", acts);
                    acts
                };
                let computed: Vec<_> = if lpt {
                    let order = x_redo_lpt_order(&wanted, |index| batch[index].2.hsps.len());
                    x_redo_lpt(pool, &states, batch.len(), &order, |worker, index| {
                        let (context, oid, list) = &batch[index];
                        redo_item(
                            *context,
                            *oid,
                            list,
                            &worker.params,
                            &mut worker.scratch,
                            &mut worker.workspace,
                            None,
                        )
                    })
                } else {
                    pool.install(|| {
                        batch
                            .par_iter()
                            .zip(wanted.par_iter())
                            .map(|((context, oid, list), &wanted)| {
                                if !wanted {
                                    return None;
                                }
                                let slot =
                                    rayon::current_thread_index().unwrap_or(0) % states.len();
                                let mut guard =
                                    states[slot].lock().expect("TBLASTN redo worker state");
                                let worker = &mut *guard;
                                Some(redo_item(
                                    *context,
                                    *oid,
                                    list,
                                    &worker.params,
                                    &mut worker.scratch,
                                    &mut worker.workspace,
                                    None,
                                ))
                            })
                            .collect()
                    })
                };
                for (at, ((context, oid, list), result)) in batch.iter().zip(computed).enumerate() {
                    let excluded = excluded_at_fetch.excluded(first + at);
                    let terminates_early = compo_early_termination(list.best_evalue, &heaps);
                    if !sw_status.reaches_redo(terminates_early, excluded) {
                        // A match that early termination stops: a redone match was done in
                        // redo_item (LOSAT_X_TBNSTAGEDROP); a marker's subject is not on the stage.
                        if terminates_early && result.is_none() && !excluded {
                            x_stage.done(*oid);
                        }
                        continue;
                    }
                    if let Some(ref mut trace) = observer {
                        trace.redo.push((
                            *context,
                            *oid,
                            list.hsps.iter().map(|hsp| hsp.hsp).collect(),
                        ));
                    }
                    let item = match result {
                        Some(result) => result,
                        None => redo_item(
                            *context,
                            *oid,
                            list,
                            &redo_params,
                            &mut scratch,
                            &mut workspace,
                            None,
                        ),
                    };
                    let item =
                        match_outcome::<R, _>(item, &mut sw_status, redo_params.do_link_hsps)?;
                    keep_item(*context, &mut heaps, item);
                }
            }
        }
    } else {
        // One thread (also a query batch searched side by side with others, on one thread each) or
        // one match: the redo runs here, outside LOSAT_X_TBNREDOLPT's rule.
        #[cfg(feature = "parallel")]
        if x_tbn_redo_lpt() {
            crate::utils::xreach::x_reach!("LOSAT_X_TBNREDOLPT", false);
        }
        for (index, (context, oid, list)) in preliminary_lists.into_iter().enumerate() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3383-3427,3525-3736
            // Each list in stream order: the early termination test on its best e-value against
            // the batch's heaps frees the list and goes on; otherwise its query's index is read.
            let excluded = excluded_at_fetch.excluded(index);
            let terminates_early = compo_early_termination(list.best_evalue, &heaps);
            // NCBI's order (`SwBatchStatus::reaches_redo`): early termination, then the fetch, which
            // excludes a marker's subject and sets the status back to 0.
            if !sw_status.reaches_redo(terminates_early, excluded) {
                // A marker's subject is not on the stage.
                if terminates_early && !excluded {
                    x_stage.done(oid);
                }
                continue;
            }
            if let Some(ref mut trace) = observer {
                trace
                    .redo
                    .push((context, oid, list.hsps.iter().map(|hsp| hsp.hsp).collect()));
            }
            let item = redo_item(
                context,
                oid,
                &list,
                &redo_params,
                &mut scratch,
                &mut workspace,
                match_redo
                    .reads_heaps(&redo_params)
                    .then_some(&mut heaps[..]),
            );
            let item = match_outcome::<R, _>(item, &mut sw_status, redo_params.do_link_hsps)?;
            keep_item(context, &mut heaps, item);
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779 (cited at
    // `match_redo::sw_batch_failed`)
    // The status of the last match the Smith-Waterman redo reached fails the batch.
    if R::SMITH_WATERMAN && sw_status.failed() {
        return Err(sw_batch_failed(
            heaps.first().is_some_and(|heap| heap.len() > 0),
        ));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2494-2515;
    // c++/src/algo/blast/core/blast_hits.c:3243-3297,3383-3400,3420-3437:
    // Blast_HSPResultsSortByEvalue(results); /* positive local TotLen */
    let mut results = Vec::with_capacity(queries.len());
    for (context, heap) in heaps.iter_mut().enumerate() {
        let mut hitlist = KappaResultHitList::new(profile.max_target_seqs)?;
        while let Some(popped) = heap.pop() {
            let list = owned_lists
                .remove(&(context, popped.subject_index))
                .context("TBLASTN composition heap lost its HSP payload")?;
            hitlist.update(list)?;
        }
        hitlist.reverse_order();
        if total_nt > 0 && profile.subject_set.sorts_results() {
            hitlist.sort_by_evalue_for_positive_totlen();
        }
        results.push(hitlist);
    }
    Ok(results)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:153-154
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:338-344
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:496-501,597
// With two or more search threads NCBI's database reader cuts the subjects into chunks by
// residues (x_FillSeqBuffer, below), not into the hundredth of the subjects of its one-thread
// iterator (`blast_engine.c:1409`). The budget of a chunk is the slice size (the largest file
// the reader opened, at most 1 GiB) over four times the threads, plus one. `-subject` has no
// database file; its F is the size of the sequence file `makeblastdb` writes for the same set,
// a quarter byte per base. LOSAT's tblastn has no `-db` yet (`args.rs` rejects it); it would
// pass the largest file of the database. Scheduling only: the number enters no computation.
/// Graduated 2026-10-22 (was `LOSAT_X_TBNOIDCHUNK`): the bases one chunk of subjects may hold.
#[cfg(feature = "parallel")]
fn x_oid_chunk_bases(threads: usize, largest_file: u64) -> u64 {
    // NCBI's Int8 `GetSliceSize() / (4*m_NumThreads) + 1`; usize fits u64 on every target.
    let divisor = u64::try_from(threads.max(1))
        .unwrap_or(u64::MAX)
        .saturating_mul(4);
    largest_file.min(1 << 30) / divisor + 1
}

// No NCBI counterpart: the size of the `.nsq` file that `makeblastdb` writes for the subjects,
// as F of `x_oid_chunk_bases` for `-subject`; it does not change any value NCBI computes.
/// ceil(bases of all subjects / 4).
#[cfg(feature = "parallel")]
fn x_subject_file_size(subjects: SubjectLetters<'_>) -> u64 {
    subjects.longest_and_total().1.div_ceil(4)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:697-707
// A chunk takes its first subject and then the next ones while the bases taken stay within
// the budget (`tot_length >= res.length` is `taken + length <= bases`). NCBI's last condition
// compares the OID after the next one with the end, so the last subject of a database always
// starts a chunk of its own, and a chunk never crosses a volume; LOSAT does not copy these, as
// `-subject` has one volume and the rest is scheduling only. NCBI reads the exact length of each
// sequence of a database from its mapped files; a `-db` stream's chunks are cut from the index's
// approximate lengths (`XDbStream::x_approx_lengths`, at most a few bases off), without a
// read of the sequence file or a table of lengths. Scheduling only.
/// The ends of the chunks of subjects, cut from their lengths given one after the other in OID
/// order (`push`).
#[cfg(feature = "parallel")]
struct XChunkCut {
    bases: u64,
    /// The bases of the chunk that the last subject joined.
    taken: u64,
    /// The subjects given.
    count: u32,
    ends: Vec<u32>,
}

#[cfg(feature = "parallel")]
impl XChunkCut {
    fn new(bases: u64) -> Self {
        Self {
            bases,
            taken: 0,
            count: 0,
            ends: Vec::new(),
        }
    }

    /// The next subject, of `length` bases: it joins the last chunk if the chunk's bases stay
    /// within the budget, else it starts the next chunk.
    fn push(&mut self, length: u64) {
        if self.count > 0 {
            match self.taken.checked_add(length) {
                Some(sum) if sum <= self.bases => self.taken = sum,
                _ => {
                    self.ends.push(self.count);
                    self.taken = length;
                }
            }
        } else {
            self.taken = length;
        }
        self.count += 1;
    }

    /// The end of every chunk, the last one's included.
    fn ends(mut self) -> Vec<u32> {
        if self.count > 0 {
            self.ends.push(self.count);
        }
        self.ends
    }
}

/// The end of the chunk of subjects that starts at `begin`: the reference of `XChunkCut`'s tests.
#[cfg(all(test, feature = "parallel"))]
fn x_oid_chunk_end(subjects: SubjectLetters<'_>, begin: usize, bases: u64) -> usize {
    let length = |oid: usize| u64::try_from(subjects.length(oid)).unwrap_or(u64::MAX);
    if begin >= subjects.count() {
        return begin;
    }
    let mut taken = length(begin);
    let mut end = begin + 1;
    while end < subjects.count() {
        match taken.checked_add(length(end)) {
            Some(sum) if sum <= bases => {
                taken = sum;
                end += 1;
            }
            _ => break,
        }
    }
    end
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:333-350
// The chunks NCBI's threads take one after another, all of them: each starts where the last one
// ended. Each NCBI search thread takes the next chunk when it is free; LOSAT's pool takes the
// chunks and each chunk is one task. Scheduling only.
/// Graduated 2026-10-22 (was `LOSAT_X_TBNOIDNEXT`, `LOSAT_X_TBNOIDCHUNK`): every chunk of the
/// subjects, in OID order (`XChunkCut`). A `-db` stream cuts its sequences once for each budget
/// and keeps the cut for the run's query batches (`XDbStream::x_kept_cut`).
#[cfg(feature = "parallel")]
fn x_all_chunks(subjects: SubjectLetters<'_>, bases: u64) -> Vec<std::ops::Range<usize>> {
    let ends: std::sync::Arc<[u32]> = match subjects {
        SubjectLetters::Stream(stream) => stream.x_kept_cut(bases, || {
            let mut cut = XChunkCut::new(bases);
            stream.x_approx_lengths(|length| cut.push(u64::from(length)));
            cut.ends()
        }),
        _ => {
            let mut cut = XChunkCut::new(bases);
            for length in subjects.lengths() {
                cut.push(u64::try_from(length).unwrap_or(u64::MAX));
            }
            cut.ends().into()
        }
    };
    let mut start = 0;
    ends.iter()
        .map(|&end| {
            let chunk = start..end as usize;
            start = end as usize;
            chunk
        })
        .collect()
}

// No NCBI counterpart: NCBI's threads write their lists into the shared HSP stream at once. The
// lists of chunks finished ahead of the oldest one LOSAT has not reduced yet wait, at most this
// many chunks beyond it, so that their memory stays bounded. Scheduling only.
/// The chunks taken or waiting beyond the oldest one not reduced.
#[cfg(feature = "parallel")]
fn x_flow_window(threads: usize) -> usize {
    threads.saturating_mul(4).max(1)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334,3346-3350,3499-3502
// One XRedoWorker is the Rust counterpart of the per-thread state NCBI allocates for the redo loop
// (NRrecord_tld, redo_align_params_tld, and the gap_align and matrix of the thread).
// One pool slot's redo state.  NCBI c++/src/algo/blast/core/blast_kappa.c:3290-3340
// allocates the same per thread (redo_align_params_tld, gap_align_tld, NRrecord_tld).
// `S` is the redo's memory (`match_redo::MatchRedo::Scratch`: NCBI's `gap_align`, or the
// Smith-Waterman redo's forbidden ranges and windows).
#[cfg(feature = "parallel")]
struct XRedoWorker<S> {
    params: BlastRedoAlignParams,
    scratch: S,
    workspace: BlastCompositionWorkspace,
}

// SAFETY: the only non-Send field is the `context` cell of the gapping
// parameters, which `redo_preliminary_match` fills for the duration of one
// synchronous call and restores before it returns; the state is only reached
// through its mutex, so it is never observed mid-call by another thread.
#[cfg(feature = "parallel")]
unsafe impl<S: Send> Send for XRedoWorker<S> {}

// No NCBI counterpart: reads the LOSAT_X_TBNREDOLPT variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNREDOLPT): a switch of the strict set (`utils/xswitch.rs`), on unless
/// it is `0`. It acts only where the redo runs on the search pool (two or more threads, two or
/// more matches), at any thread count, on a chunk with two or more matches to redo. (A limit of
/// 16 threads, set while the reap of `kappa.rs` was quadratic, came off with `LOSAT_X_TBNREAPIDX`:
/// `docs/agents/diagnostics.md`.)
#[cfg(feature = "parallel")]
fn x_tbn_redo_lpt() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNREDOLPT"))
}

// No NCBI counterpart: the order in which the threads take a chunk's matches; it changes no value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNREDOLPT): the indices of the `wanted` matches by `hsps(index)` (the
/// number of HSPs of the match's preliminary list), largest first; matches with as many HSPs in
/// index order. A redo's cost grows with its list, whatever the function, so any such cost gives
/// this order.
#[cfg(feature = "parallel")]
fn x_redo_lpt_order(wanted: &[bool], hsps: impl Fn(usize) -> usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..wanted.len()).filter(|&index| wanted[index]).collect();
    order.sort_by_key(|&index| (std::cmp::Reverse(hsps(index)), index));
    order
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3447-3449,3493-3497
// No NCBI counterpart for the order: NCBI's threads take fixed shares of its loop; here the
// threads take the matches one at a time, largest first. Scheduling only.
/// EXPERIMENT (LOSAT_X_TBNREDOLPT): `run(worker, index)` for every index of `order`, in one
/// `rayon::scope` on `pool` of `min(workers, order)` tasks, task `t` with `workers[t]`. Each task
/// takes the next index of `order` from one counter until none is left, and puts each result in
/// the index's slot. Returns the results of a chunk of `matches` matches by index, `None` where
/// `order` does not name the index (a match early termination ends before the chunk).
///
/// Why the results are the `par_iter` path's: a match's result is a function of its preliminary
/// list (`redo_item` with the same arguments); a worker carries no value from one match to the
/// next (the `par_iter` path also gives one worker the matches of whichever thread holds it); the
/// results are read by index, so the in-order apply sees the same sequence, errors included. The
/// order only changes which thread computes what, and when; no floating-point operation changes
/// order, the matches share no reduction. A task holds its own worker and waits for nothing but a
/// `-db` sequence another thread is reading (`XDbStream::read_shared`), a read that runs no rayon
/// work; a match runs none either, so no task depends on a task that could be stolen below it.
#[cfg(feature = "parallel")]
#[inline(never)]
fn x_redo_lpt<W, R>(
    pool: &SearchPool<'_>,
    workers: &[std::sync::Mutex<W>],
    matches: usize,
    order: &[usize],
    run: impl Fn(&mut W, usize) -> R + Sync,
) -> Vec<Option<R>>
where
    W: Send,
    R: Send + Sync,
{
    // A disabled pool would run the scope on the global pool, not the search pool.
    debug_assert!(pool.enabled());
    let slots: Vec<std::sync::OnceLock<R>> =
        (0..matches).map(|_| std::sync::OnceLock::new()).collect();
    let next = AtomicUsize::new(0);
    let (slots_ref, next, run) = (&slots, &next, &run);
    pool.install(|| {
        rayon::scope(|scope| {
            for worker in &workers[..workers.len().min(order.len())] {
                scope.spawn(move |_| {
                    let mut guard = worker.lock().expect("TBLASTN redo worker state");
                    while let Some(&index) = order.get(next.fetch_add(1, Ordering::Relaxed)) {
                        let result = run(&mut guard, index);
                        assert!(
                            slots_ref[index].set(result).is_ok(),
                            "match {index} ran twice"
                        );
                    }
                });
            }
        });
    });
    slots
        .into_iter()
        .map(std::sync::OnceLock::into_inner)
        .collect()
}

// The per-context cutoffs of the preliminary stage of a query set.
struct ContextCutoffArrays {
    word_xdrop: Vec<i32>,
    word_cutoff: Vec<i32>,
    hit_cutoff: Vec<i32>,
}

impl ContextCutoffArrays {
    fn of(parameters: &LocalSubjectParameters) -> Self {
        Self {
            word_xdrop: parameters.cutoffs.iter().map(|v| v.word_xdrop).collect(),
            word_cutoff: parameters.cutoffs.iter().map(|v| v.word_cutoff).collect(),
            hit_cutoff: parameters.cutoffs.iter().map(|v| v.hit_cutoff).collect(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:245 and 1292
/// The word threshold in the units of the lookup table: the standard table's (word sizes up to
/// 4), or the compressed table's, scaled by `kMatrixScale` 100 (word sizes 5 to 7, the only ones
/// that take it: `SetWordSize`, blast_options_local_priv.hpp:612-619), converted as the x86-64
/// build does (`INT_MIN` for a value beyond Int4 or not a number).
pub(super) fn lookup_threshold(word_size: usize, threshold: f64) -> i32 {
    let scale = if word_size > 4 { 100.0 } else { 1.0 };
    crate::core::blast_util::ncbi_int4_from_double(scale * threshold)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:302-383,455-463:
// The word and extension parameters of the preliminary stage of a query set.
fn preliminary_profile<'a>(
    profile: LocalStageDProfile<'a>,
    scoring: LocalStageDScoring,
    gapped: KarlinParams,
    cutoffs: &'a ContextCutoffArrays,
) -> PreliminaryProfile<'a> {
    PreliminaryProfile {
        seg: profile.seg,
        soft_masking: profile.soft_masking,
        threshold: lookup_threshold(scoring.word_size, scoring.threshold),
        window: scoring.window,
        word_xdrop: &cutoffs.word_xdrop,
        word_cutoff: &cutoffs.word_cutoff,
        mask_lowercase: profile.mask_lowercase,
        matrix: scoring.matrix,
        word_size: scoring.word_size,
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:457-458: (Int4)(gap_x_dropoff*NCBIMATH_LN2/min_lambda),
        // converted as the x86-64 build does (INT_MIN beyond Int4; `core::x_drop`).
        gap_xdrop: crate::core::x_drop::gapped_x_dropoffs(
            scoring.gap_xdrop_bits,
            scoring.final_xdrop_bits,
            gapped.lambda,
        )
        .0,
        gapped_cutoff: &cutoffs.hit_cutoff,
        hsp_num_max: i32::MAX as usize,
        subject_besthit: profile.subject_besthit,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-290
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:171-182
// The preliminary stage of a split batch (RP-4): each chunk is a query set of its parts
// (the chunk's query factory, with the lower-case masks restricted to each part), with the
// batch's effective search spaces, and its hit lists merge into the batch's. A chunk whose
// contexts all fail the Karlin-Altschul setup finds nothing (NCBI ignores its error).
#[allow(clippy::too_many_arguments)]
fn split_preliminary_hitlists(
    queries: &[Vec<u8>],
    chunks: &[ProteinQueryChunk],
    sizes: SplitSizes,
    batch_parameters: &LocalSubjectParameters,
    option_spaces: Option<&[i64]>,
    inputs: QuerySetInputs<'_>,
    subjects: SubjectLetters<'_>,
    prelim_size: usize,
    pool: Option<&SearchPool<'_>>,
) -> Result<Vec<Option<KappaResultHitList>>> {
    // A negative CHUNK_SIZE splits a batch only with a negative OVERLAP_CHUNK_SIZE below
    // it; NCBI's `size_t` chunk ranges then wrap (split_query_cxx.cpp:145-171).
    ensure!(
        !sizes.negative_chunk_size(),
        "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this query batch with NCBI BLAST+'s size_t chunk ranges wrapped, which is not supported by LOSAT's TBLASTN"
    );
    ensure!(
        !protein_chunk_would_be_split(chunks, sizes),
        "an OVERLAP_CHUNK_SIZE so close to CHUNK_SIZE that NCBI BLAST+ would split a query chunk again (it then stops with a null reference) is not supported by LOSAT's TBLASTN"
    );
    // The options' spaces of the chunks (`SplitQuery_SetEffectiveSearchSpace`): the batch's own,
    // or the preliminary stage's where they differ (`run_local_search_with_pool`).
    let batch_eff_searchsp: Vec<i64> = match option_spaces {
        Some(spaces) => spaces.to_vec(),
        None => batch_parameters
            .lengths
            .iter()
            .map(|length| length.eff_searchsp)
            .collect(),
    };
    let per_subject = inputs.profile.subject_set.per_subject();
    let mut merged: Vec<Option<KappaResultHitList>> = (0..queries.len()).map(|_| None).collect();
    for chunk in chunks {
        // NCBI fails to make the query factory of a chunk without a query ("Empty
        // CBlastQueryVector"); `split_protein_batch` makes none.
        ensure!(
            !chunk.queries.is_empty(),
            "a query chunk without a query is not supported by LOSAT's TBLASTN"
        );
        let parts: Vec<Vec<u8>> = chunk
            .queries
            .iter()
            .map(|part| {
                let query = &queries[part.query];
                if !inputs.profile.mask_lowercase {
                    return query[part.from..part.to].to_vec();
                }
                let (record, offset) = inputs
                    .profile
                    .query_records
                    .get(part.query)
                    .copied()
                    .unwrap_or((query, 0));
                protein_chunk_part_in_record(record, offset, part)
            })
            .collect();
        // A chunk's setup is inside `search` (the common line marks only the batch's own).
        let setup = query_set_setup(&parts, inputs, Some(&batch_eff_searchsp), false, pool)?;
        if !setup.valid_contexts.iter().any(|&valid| valid) {
            continue;
        }
        let query_lengths: Vec<i32> = parts
            .iter()
            .map(|part| i32::try_from(part.len()))
            .collect::<std::result::Result<_, _>>()?;
        let cutoff_arrays = ContextCutoffArrays::of(&setup.parameters);
        let part_refs: Vec<&[u8]> = parts.iter().map(Vec::as_slice).collect();
        let chunk_lists = run_preliminary_subjects(
            &PreliminaryInputs {
                queries: &part_refs,
                frames: &setup.working_frames,
                genetic_code: inputs.profile.genetic_code,
                search: preliminary_profile(
                    inputs.profile,
                    inputs.scoring,
                    inputs.gapped,
                    &cutoff_arrays,
                ),
                link: setup.parameters.link,
                query_lengths: &query_lengths,
                lengths: &setup.parameters.lengths,
                gapped_by_context: &setup.gapped_by_context,
                gumbel: inputs.gumbel,
                prelim_evalue: setup.parameters.prelim_evalue,
                prelim_size,
                culling_limit: inputs.profile.culling_limit,
                pssm: inputs.profile.pssm.map(|psi| &psi.pssm),
                update: per_subject.then(|| SubjectUpdate {
                    recipe: &setup.recipe,
                    option_spaces: &batch_eff_searchsp,
                    dbsize_nt: inputs.profile.subject_set.dbsize_nt(),
                    min_length: inputs.scoring.word_size * 3 + 2,
                }),
            },
            subjects,
            pool,
            None,
        )?;
        let parts_and_offsets: Vec<(usize, i32)> = chunk
            .queries
            .iter()
            .zip(&chunk.context_offsets)
            .map(|(part, &offset)| (part.query, offset))
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:174-178
        merge_query_chunk(
            &mut merged,
            chunk_lists,
            &parts_and_offsets,
            sizes.overlap as u32 as i32,
        )?;
    }
    Ok(merged)
}

// The query-set inputs of `query_set_setup`: the whole query batch, or the parts of one query
// chunk of a split batch (RP-4).
#[derive(Clone, Copy)]
struct QuerySetInputs<'a> {
    profile: LocalStageDProfile<'a>,
    scoring: LocalStageDScoring,
    gapped: KarlinParams,
    ungapped: &'a KarlinParams,
    gumbel: Option<&'a BlastGumbelBlk>,
    /// The database length of the effective lengths in nucleotides
    /// (`SubjectSetNumbers::eff_length_nt`), NCBI's `Int8`.
    total_nt: i64,
    num_subjects: i32,
    /// The setup's subject length of the hit-saving parameters (`setup_subject_length`), 0
    /// where the setup fails (`SetupEngineError::SUBJECT_LENGTH_INVALID`).
    min_subject_length: i32,
    compo_adjust_mode: BlastCompoAdjustMode,
    do_sum_stats: bool,
}

struct QuerySetSetup<'a> {
    working_frames: Vec<Vec<crate::algorithm::tblastx::translation::QueryFrame>>,
    contexts: Vec<crate::algorithm::tblastx::lookup::QueryContext>,
    valid_contexts: Vec<bool>,
    gapped_by_context: Vec<KarlinParams>,
    parameters: LocalSubjectParameters,
    /// What the parameters of the query set are computed from, for the preliminary stage's own
    /// (`LocalParameterRecipe::parameters`).
    recipe: LocalParameterRecipe<'a>,
}

/// The inputs of `local_parameters_for_blocks` of a query set (`query_set_setup`): its contexts,
/// their Karlin-Altschul blocks and the options, with which the setup's parameters, the preliminary
/// stage's (its option spaces: `stats::setup::preliminary_option_space`) and those of each subject
/// (`LocalParameterCall::Update`) are computed.
struct LocalParameterRecipe<'a> {
    query_contexts: Vec<(usize, bool)>,
    gapped: Vec<KarlinParams>,
    lengths: Vec<KarlinParams>,
    ungapped: Vec<KarlinParams>,
    x_dropoff: Vec<KarlinParams>,
    num_subjects: i32,
    options: LocalParameterOptions<'a>,
}

impl LocalParameterRecipe<'_> {
    /// The parameters with the database length `db_nt` in nucleotides (the setup's total, `-dbsize`
    /// or a subject's length; NCBI's `Int8` `real_db_length`), the options' spaces `option_spaces`
    /// and the call `call`.
    fn parameters(
        &self,
        db_nt: i64,
        option_spaces: Option<&[i64]>,
        call: LocalParameterCall,
    ) -> Result<LocalSubjectParameters, crate::stats::setup::SearchSpaceOverflow> {
        super::stage_d_stats::local_parameters_for_blocks(
            &self.query_contexts,
            db_nt,
            self.num_subjects,
            super::stage_d_stats::LocalKarlinBlocks {
                gapped: &self.gapped,
                lengths: &self.lengths,
                ungapped: &self.ungapped,
                x_dropoff: &self.x_dropoff,
            },
            LocalParameterOptions {
                eff_searchsp_override: option_spaces,
                ..self.options
            },
            call,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-654,964-1024:
// BlastHitSavingParametersNew(...); the score block, effective lengths and cutoffs of a
// query set. `eff_searchsp_override` is the options' effective search spaces that NCBI sets
// from the whole batch before it searches a query chunk (`SplitQuery_SetEffectiveSearchSpace`).
fn query_set_setup<'a>(
    queries: &[Vec<u8>],
    inputs: QuerySetInputs<'a>,
    eff_searchsp_override: Option<&[i64]>,
    marked: bool,
    pool: Option<&SearchPool<'_>>,
) -> Result<QuerySetSetup<'a>> {
    // No NCBI counterpart (scheduling only): NCBI masks the queries of a batch and computes their
    // contexts' Karlin blocks one after the other on the thread that sets the batch up. Each reads
    // its own query alone, so a batch searched with a pool of two or more threads inside it (one
    // after another, `run_local_for_report_threads`) makes them on the pool's threads, where they
    // would otherwise hold the others idle; a batch searched side by side runs them on its thread.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let x_pool = pool.filter(|pool| pool.enabled() && pool.threads() > 1);
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let _ = pool;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2803;
    // c++/src/algo/blast/core/blast_setup.c:770-847:
    // Reuse the Stage C score-frequency calculation on the query bytes after
    // BlastSetUp_MaskQuery, which precedes BlastSetup_ScoreBlkInit at
    // pinned blast_setup.c:614-654.
    let encode = |query: &Vec<u8>, timed: bool| {
        // Diagnostic (LOSAT_TIMING): with `timed`, the query's encoding and SEG go to the
        // common line (utils/xtiming.rs). No NCBI counterpart.
        let frame = if inputs.profile.soft_masking {
            if timed {
                encode_protein_query_frame_with_seg_timed(query, None)
            } else {
                encode_protein_query_frame_with_seg(query, None)
            }
        } else if timed {
            super::search_seed::encode_tblastn_lookup_query_timed(
                query,
                inputs.profile.seg,
                inputs.profile.mask_lowercase,
            )
        } else {
            super::search_seed::encode_tblastn_lookup_query(
                query,
                inputs.profile.seg,
                inputs.profile.mask_lowercase,
            )
        };
        vec![frame]
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let working_frames: Vec<_> = match x_pool {
        Some(pool) => {
            use rayon::prelude::*;
            // Diagnostic (LOSAT_TIMING): with `marked`, the queries' encoding and SEG on the pool
            // as one interval of this thread (`seg`; the timed encoders mark the calling thread's
            // own intervals). No NCBI counterpart.
            let t_seg = if marked {
                crate::utils::xtiming::clock()
            } else {
                None
            };
            let frames = pool.install(|| queries.par_iter().map(|q| encode(q, false)).collect());
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
            frames
        }
        None => queries.iter().map(|q| encode(q, marked)).collect(),
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let working_frames: Vec<_> = queries.iter().map(|q| encode(q, marked)).collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2803;
    // c++/src/algo/blast/core/lookup_wrap.c:91-100:
    // The active matrix determines context validity; the active word size
    // determines the lookup state used by Stage C.
    // Only the contexts are used here (graduated 2026-10-22, was LOSAT_X_LUTSPLIT).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2792
    // Only the contexts (their validity, from the ungapped Karlin block) are used here, not the
    // lookup table. x_lookup_contexts_for_profile returns the contexts without building the table,
    // which the search builds where it needs it (lookup_wrap.c:91-100).
    // Diagnostic (LOSAT_TIMING): the query's lookup contexts for the common line
    // (utils/xtiming.rs), with `marked`. The table the scan reads is built inside the search
    // (before its subjects on a pool of two or more threads, `XQuerySide::x_build_lookup_on`,
    // else by the first subject that needs it, `x_scan_plan`), so it is part of `search`. No NCBI
    // counterpart.
    let t_build_lookup = if marked {
        crate::utils::xtiming::clock()
    } else {
        None
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let x_contexts = match x_pool {
        Some(pool) => crate::algorithm::tblastx::lookup::x_lookup_contexts_for_profile_in_pool(
            &working_frames,
            false,
            inputs.scoring.matrix,
            inputs.scoring.word_size,
            pool,
        ),
        None => crate::algorithm::tblastx::lookup::x_lookup_contexts_for_profile(
            &working_frames,
            false,
            inputs.scoring.matrix,
            inputs.scoring.word_size,
        ),
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let x_contexts = crate::algorithm::tblastx::lookup::x_lookup_contexts_for_profile(
        &working_frames,
        false,
        inputs.scoring.matrix,
        inputs.scoring.word_size,
    );
    let mut contexts = match x_contexts {
        Some(contexts) => contexts,
        // The compressed table (word sizes 5 to 7) is not the standard table that the diagnostics
        // count; the contexts are the same.
        None if inputs.scoring.word_size > 4 => {
            crate::algorithm::tblastx::lookup::x_compressed_contexts_for_profile(
                &working_frames,
                false,
                inputs.scoring.matrix,
                inputs.scoring.word_size,
            )
        }
        None => {
            build_ncbi_lookup_for_profile(
                &working_frames,
                lookup_threshold(inputs.scoring.word_size, inputs.scoring.threshold),
                inputs.ungapped,
                false,
                inputs.scoring.matrix,
                inputs.scoring.word_size,
            )
            .1
        }
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::BuildLookup, t_build_lookup);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2772
    // A query without letters (`Sequence contains no data`) is invalidated when the queries
    // are set up, before its Karlin-Altschul parameters would be computed: it is never
    // searched, and its effective lengths stay 0 (`local_parameters_for_call`).
    for context in &mut contexts {
        if context.aa_len == 0 {
            context.is_valid = false;
            context.karlin_params = Default::default();
        }
    }
    let query_contexts: Vec<_> = queries
        .iter()
        .zip(&contexts)
        .map(|(query, context)| (query.len(), context.is_valid))
        .collect();
    let valid_contexts: Vec<_> = contexts.iter().map(|context| context.is_valid).collect();
    let gapped_by_context = vec![inputs.gapped; queries.len()];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2777-2827;
    // c++/src/algo/blast/core/blast_parameters.c:208-223,335-347:
    // Word x-drop and gap trigger both use each query context's Karlin block.
    let ungapped_by_context: Vec<_> = contexts
        .iter()
        .map(|context| context.karlin_params)
        .collect();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:964-985,1011-1024;
    // c++/src/algo/blast/core/blast_parameters.c:457-463:
    let options = LocalParameterOptions {
        expect_value: inputs.profile.expect_value,
        do_sum_stats: inputs.do_sum_stats,
        max_intron_length: inputs.profile.max_intron_length,
        gap_trigger_bits: 22.0,
        // `-xdrop_ungap` or 7 bits (`LocalStageDScoring::word_xdrop_bits`).
        word_xdrop_bits: inputs.scoring.word_xdrop_bits,
        scale_factor: 1.0,
        gumbel: inputs.gumbel,
        eff_searchsp_override: None,
    };
    let call = LocalParameterCall::Initial {
        min_subject_length: inputs.min_subject_length,
        composition_based_stats: inputs.compo_adjust_mode.ncbi_value(),
    };
    let recipe = match inputs.profile.pssm {
        None => LocalParameterRecipe {
            query_contexts,
            gapped: gapped_by_context.clone(),
            lengths: gapped_by_context,
            ungapped: ungapped_by_context.clone(),
            x_dropoff: ungapped_by_context,
            num_subjects: inputs.num_subjects,
            options,
        },
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:89-127 (cited at
        // `PsiSearch::ungapped_karlin` and `PsiSearch::gapped_karlin`)
        // A PSSM query's `kbp_psi` and `kbp_gap_psi`: the file's values, else `kbp_std`'s and
        // `kbp_gap_std`'s. They are `kbp` and `kbp_gap` from here on (the word X-drop, the
        // cutoffs, the linking, the X-drops of the gapped stages, the e-values and bit scores);
        // the effective lengths keep `kbp_gap_std` and the gap trigger `kbp_std`
        // (`LocalKarlinBlocks`).
        Some(psi) => {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:960-961
            // `BlastExtensionParametersNew` looks for a valid `sbp->kbp`, which is `kbp_psi`
            // (blast_stat.c:2826; cited at `SetupEngineError::NO_VALID_KARLIN_BLOCK`), before the
            // average subject length (`SUBJECT_LENGTH_INVALID`, below). A PSSM query whose letters give its
            // context a valid `kbp_std` but whose file's ungapped values make `kbp_psi` invalid
            // (`lambdaUngapped`, `kappaUngapped` or `hUngapped` negative) stops the search's
            // setup, before any subject is read, with NCBI's engine error of an invalid query
            // (probes `x.ungneg.*`, `r3.*ungneg.*` of the PS7 review). A PSSM query without a
            // valid context is the one that NCBI crashes with, refused after the batches.
            if valid_contexts.iter().any(|&valid| valid)
                && !contexts.iter().any(|context| {
                    context.is_valid && psi.ungapped_karlin_valid(&context.karlin_params)
                })
            {
                return Err(SetupEngineError::NO_VALID_KARLIN_BLOCK.into());
            }
            LocalParameterRecipe {
                query_contexts,
                gapped: gapped_by_context
                    .iter()
                    .map(|kbp| psi.gapped_karlin(kbp))
                    .collect(),
                lengths: gapped_by_context,
                x_dropoff: ungapped_by_context
                    .iter()
                    .map(|kbp| psi.ungapped_karlin(kbp))
                    .collect(),
                ungapped: ungapped_by_context,
                num_subjects: inputs.num_subjects,
                options,
            }
        }
    };
    let parameters = recipe.parameters(inputs.total_nt, eff_searchsp_override, call)?;
    let gapped_by_context = recipe.gapped.clone();
    Ok(QuerySetSetup {
        working_frames,
        contexts,
        valid_contexts,
        gapped_by_context,
        parameters,
        recipe,
    })
}

// The inputs of the preliminary stage over the subjects for one query set (the whole batch,
// or one query chunk of a split batch).
struct PreliminaryInputs<'a> {
    queries: &'a [&'a [u8]],
    /// The queries' frames, masked by the setup of the query set (`query_set_setup`).
    frames: &'a [Vec<crate::algorithm::tblastx::translation::QueryFrame>],
    genetic_code: u8,
    search: PreliminaryProfile<'a>,
    link: Option<super::stage_d_stats::LocalLinkParameters>,
    query_lengths: &'a [i32],
    lengths: &'a [super::stage_d_stats::LocalContextLength],
    gapped_by_context: &'a [KarlinParams],
    gumbel: Option<&'a BlastGumbelBlk>,
    prelim_evalue: f64,
    prelim_size: usize,
    /// `-culling_limit` (0: none).
    culling_limit: i32,
    /// `-in_pssm`: the PSSM's rows, the scores of the one query (`score::PssmScore`) and the
    /// source of its lookup table (`in_pssm::lookup_table`); `None` for a FASTA batch.
    pssm: Option<&'a crate::algorithm::pssm_search::rows::QueryPssm>,
    /// The statistics of each subject (`SubjectSetNumbers::per_subject`), `None` for the query
    /// set's.
    update: Option<SubjectUpdate<'a>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1429-1443
/// What the statistics of each subject read (`PreliminaryInputs::update`): the query set's
/// parameter inputs, the preliminary stage's option spaces of its contexts, `-dbsize` in
/// nucleotides when above 0 (the effective lengths' database length in place of the subject's),
/// and NCBI's shortest searched subject (`min_subj_seq_length`, word size * 3 + 2), which a shorter
/// subject does not reach.
#[derive(Clone, Copy)]
struct SubjectUpdate<'a> {
    recipe: &'a LocalParameterRecipe<'a>,
    option_spaces: &'a [i64],
    dbsize_nt: Option<i64>,
    min_length: usize,
}

// Where the HSP lists of the subjects go as the owner of the search takes them in OID order:
// the hit-list collector of the preliminary stage (`PreliminaryCollector`), or, with
// `-culling_limit`, NCBI's culling writer (`PrelimCulling`). `run_preliminary_subjects` picks
// one once per call, before the subject loop. The search takes the sink as a trait object, so
// it has one instance (one indirect call per subject), as before the culling writer.
trait PreliminarySink {
    /// The list of subject `oid`: the linked HSPs of all queries and frames, in list order.
    fn add(&mut self, oid: usize, linked: &LinkedHspList) -> Result<()>;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:129-161:
// Each context owns an independently capped preliminary list.
struct PreliminaryCollector {
    hitlists: Vec<Option<KappaResultHitList>>,
    prelim_size: usize,
}

impl PreliminaryCollector {
    fn new(inputs: &PreliminaryInputs<'_>) -> Self {
        Self {
            hitlists: (0..inputs.queries.len()).map(|_| None).collect(),
            prelim_size: inputs.prelim_size,
        }
    }
}

impl PreliminarySink for PreliminaryCollector {
    fn add(&mut self, oid: usize, linked: &LinkedHspList) -> Result<()> {
        for context in 0..self.hitlists.len() {
            let hsps: Vec<_> = linked
                .hsps
                .iter()
                .copied()
                .filter(|h| h.context == context)
                .collect();
            if !hsps.is_empty() {
                let best_evalue = hsps.iter().map(|h| h.evalue).reduce(f64::min).unwrap();
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:129-161:
                // Each context owns an independently capped preliminary list.
                let collector = self.hitlists[context]
                    .get_or_insert(KappaResultHitList::new(self.prelim_size)?);
                collector.update(KappaResultList {
                    oid: i32::try_from(oid)?,
                    hsps: LinkedHspList { hsps, best_evalue },
                    payloads: Vec::new(),
                })?;
            }
        }
        Ok(())
    }
}

impl PreliminarySink for PrelimCulling {
    fn add(&mut self, oid: usize, linked: &LinkedHspList) -> Result<()> {
        PrelimCulling::add(self, oid, &linked.hsps)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1414,1469-1475,1553-1555:
// The preliminary hit list of each query context of the query set.
fn run_preliminary_subjects(
    inputs: &PreliminaryInputs<'_>,
    subjects: SubjectLetters<'_>,
    pool: Option<&SearchPool<'_>>,
    observer: Option<&mut StageDBoundaryTrace>,
) -> Result<Vec<Option<KappaResultHitList>>> {
    // `-culling_limit` selects NCBI's culling writer in place of the collector, once per call.
    if inputs.culling_limit > 0 {
        let mut sink = PrelimCulling::new(
            inputs.query_lengths,
            inputs.prelim_size,
            inputs.culling_limit,
        );
        run_preliminary_subjects_into(inputs, subjects, pool, observer, &mut sink)?;
        sink.finish()
    } else {
        let mut sink = PreliminaryCollector::new(inputs);
        run_preliminary_subjects_into(inputs, subjects, pool, observer, &mut sink)?;
        Ok(sink.hitlists)
    }
}

fn run_preliminary_subjects_into(
    inputs: &PreliminaryInputs<'_>,
    subjects: SubjectLetters<'_>,
    pool: Option<&SearchPool<'_>>,
    mut observer: Option<&mut StageDBoundaryTrace>,
    sink: &mut (dyn PreliminarySink + Send),
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905,1469-1475:
    // Finish search, link/direct E-value and preliminary reap inside one subject
    // job, before its OID-ordered list enters the shared collector.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
    // NCBI scans the subjects in this loop. With one job per subject the pool is busy with
    // subjects; otherwise the units of the single subject use it. Scheduling only.
    // Each subject keeps the NCBI search-to-link-to-reap order. Indexed
    // collection finishes before the OID-ordered shared collector below.
    let parallel_selected = pool.is_some_and(SearchPool::enabled) && subjects.count() > 1;
    let unit_pool = if parallel_selected { None } else { pool };
    // Graduated 2026-10-22 (was `LOSAT_X_TBNOIDCHUNK`): the bases one chunk may hold, from the
    // size of the sequence file of the subjects (x_oid_chunk_bases); none when the subjects are
    // not searched on the pool.
    #[cfg(feature = "parallel")]
    crate::utils::xreach::x_reach!("LOSAT_X_TBNOIDCHUNK", parallel_selected);
    #[cfg(feature = "parallel")]
    let x_chunk_bases = match pool {
        Some(pool) if parallel_selected => Some(x_oid_chunk_bases(
            pool.threads(),
            x_subject_file_size(subjects),
        )),
        _ => None,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-482
    // Graduated 2026-10-22 (was `LOSAT_X_TBNOIDUNITS`): NCBI's thread searches the chunks of its
    // subject one after another. A subject longer than a chunk's budget, a chunk of its own, runs
    // its (frame, chunk) units on the pool as the single subject of a search does, so that threads
    // without a chunk to search take them; its merge, linking and reaping stay on its thread after
    // its last unit. The other subjects keep `unit_pool`. Scheduling only. (A `-db` stream's chunks
    // are cut from the index's lengths, a few bases off the subject's length that decides here:
    // a subject within those bases of the budget may share its chunk and still run its units on
    // the pool.)
    #[cfg(feature = "parallel")]
    let x_units_beyond = x_chunk_bases;
    let x_unit_pool = |subject: XPrelimSubject<'_>| {
        #[cfg(feature = "parallel")]
        if let Some(bases) = x_units_beyond {
            let beyond = u64::try_from(subject.len()).unwrap_or(u64::MAX) > bases;
            crate::utils::xreach::x_reach!("LOSAT_X_TBNOIDUNITS", beyond);
            if beyond {
                return pool;
            }
        }
        #[cfg(not(feature = "parallel"))]
        let _ = subject;
        unit_pool
    };
    // The lookup table of this query set, built by the first subject that needs it, from the
    // query frames the setup masked (`query_set_setup`; NCBI masks the batch's query once, in
    // BLAST_MainSetUp).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100
    // NCBI builds the lookup table once per query batch; this value lets the subjects of the batch
    // share the query frames and the table (`XQuerySide`).
    let x_query_side =
        super::search_gapped::XQuerySide::x_with_frames(inputs.frames, inputs.queries)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:86-99 (cited at
    // `pssm_search::lookup`)
    // A PSSM query's table, of its rows over the complement of its masks, made with the batch
    // (`in_pssm::lookup_table`).
    if let (Some(rows), [query]) = (inputs.pssm, inputs.queries) {
        x_query_side.x_set_lookup(std::sync::Arc::new(super::in_pssm::lookup_table(
            rows,
            query,
            inputs.search.seg,
            inputs.search.mask_lowercase,
            inputs.search.word_size,
            inputs.search.threshold,
        )));
    }
    // A batch searched with a pool of two or more threads inside it builds its table before its
    // subjects, on the pool (`XQuerySide::x_build_lookup_on`): the subjects would otherwise wait
    // while the first one builds it on its thread. A batch on one thread (side by side, or
    // `-num_threads 1`) builds it at its first subject that needs it. Scheduling only.
    if let Some(pool) = pool.filter(|pool| pool.enabled() && pool.threads() > 1) {
        x_query_side.x_build_lookup_on(inputs.search, pool);
    }
    // The subject plans (was LOSAT_X_TBNSSIDE) the owner of the search
    // attached to the pool; `subjects` is its subject list.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
    // The subject plans (translations) attached to the pool by the owner of the search
    // (`XSubjectSide`); a reuse of the translation NCBI makes for each subject.
    let x_subject_side = pool.and_then(|pool| pool.x_ext::<super::search_seed::XSubjectSide>());
    // No NCBI counterpart: x_oid is only the index of the subject among the shared subject plans;
    // it does not change any value NCBI computes.
    let subject_core = |x_oid: usize,
                        subject: XSearchSubject<'_>|
     -> Result<(Vec<(usize, GappedHsp)>, LinkedHspList)> {
        // A `-db` sequence's length (the linking's and the E-values' subject length below, and
        // the plan's): the one an earlier batch gave its kept plan (`XSubjectSide::length`, not
        // read again), else read where the sequence is read (`XDbStream::length`: a short one's
        // block, which its plan then reads; a long one's index and last byte). A read that fails
        // fails the subject as a read of its plan does.
        let subject = match subject {
            XSearchSubject::Nt(subject) => XPrelimSubject::Nt(subject),
            XSearchSubject::Stream(stream, oid) => {
                let kept = x_subject_side.and_then(|side| {
                    side.length(oid, inputs.genetic_code, inputs.search.mask_lowercase)
                });
                let length = match kept {
                    Some(length) => length,
                    None => stream.length(oid)?,
                };
                XPrelimSubject::Stream(stream, oid, length)
            }
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1412-1413,1469-1470
        // Dispatch point: the call that runs the preliminary stage for one subject
        // (s_BlastSearchEngineCore). preliminary_protein_hsps_for_search runs the unit-wise stage.
        let x_pool = x_unit_pool(subject);
        // The statistics of each subject that NCBI searches (`SubjectUpdate`): the query set's
        // parameters with the subject's length in nucleotides (`LocalParameterCall::Update`).
        let x_update;
        let x_arrays;
        let x_inputs;
        let inputs = match inputs.update {
            Some(update) if subject.len() >= update.min_length => {
                x_update = update.recipe.parameters(
                    update.dbsize_nt.unwrap_or(subject.len() as i64),
                    Some(update.option_spaces),
                    LocalParameterCall::Update {
                        subject_length: i32::try_from(subject.len())?,
                    },
                )?;
                x_arrays = ContextCutoffArrays::of(&x_update);
                x_inputs = PreliminaryInputs {
                    search: PreliminaryProfile {
                        word_cutoff: &x_arrays.word_cutoff,
                        gapped_cutoff: &x_arrays.hit_cutoff,
                        ..inputs.search
                    },
                    link: x_update.link,
                    lengths: &x_update.lengths,
                    prelim_evalue: x_update.prelim_evalue,
                    ..*inputs
                };
                &x_inputs
            }
            _ => inputs,
        };
        // Dispatch point: the scores of the batch's PSSM query (`-in_pssm`, the same for every
        // subject), else the matrix's.
        let preliminary = match inputs.pssm {
            None => super::search_gapped::preliminary_protein_hsps_for_search(
                inputs.queries,
                subject,
                inputs.genetic_code,
                inputs.search,
                x_pool,
                &x_query_side,
                x_subject_side.map(|side| (side, x_oid)),
                super::score::MatrixScore,
            )?,
            Some(rows) => super::search_gapped::preliminary_protein_hsps_for_search(
                inputs.queries,
                subject,
                inputs.genetic_code,
                inputs.search,
                x_pool,
                &x_query_side,
                x_subject_side.map(|side| (side, x_oid)),
                super::score::PssmScore(rows),
            )?,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905:
        // Even an allocated empty subject list enters link/E-value. The query contexts of a
        // subject whose units ran on the pool are linked there side by side (LOSAT_X_TBNLINKCTX).
        let mut linked = if let Some(link) = inputs.link {
            link_preliminary_hsps_on(
                x_pool,
                &preliminary,
                inputs.query_lengths,
                inputs.lengths,
                i32::try_from(subject.len())?,
                inputs.gapped_by_context,
                inputs.gumbel,
                &link,
            )?
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:728-729,804-813,882-887:
            let stat_length = initial_translated_stat_length(subject.len());
            direct_local_evalues(
                &preliminary,
                inputs.query_lengths,
                i32::try_from(stat_length)?,
                inputs.gapped_by_context,
                inputs.lengths,
                inputs.gumbel,
            )
        };
        reap_by_evalue(&mut linked, inputs.prelim_evalue);
        Ok((preliminary, linked))
    };
    // Subject `x_oid` for its preliminary stage: the record's letters, or (EXPERIMENT,
    // LOSAT_X_TBNSTORE) its packed form, translated from its packed bases as NCBI translates a
    // subject from ncbi2na (the letters made for the time of the call if the stage reads the case).
    let subject_core_oid = |x_oid: usize| {
        subjects.with_search_subject(x_oid, inputs.search.mask_lowercase, |subject| {
            subject_core(x_oid, subject)
        })
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1414,1469-1475:
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-188:
    crate::utils::threading::report_stage(
        "tblastn",
        "subject_core",
        subjects.count(),
        parallel_selected,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-145:
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1414:
    // Keep at most one pool-sized batch of subject results before OID reduction;
    // memory follows NCBI's bounded number of workers, not total subject count.
    // The subjects are cut into chunks by residues as NCBI's database reader cuts them for its
    // search threads (x_all_chunks), and the flow of chunks (`x_oid_flow.rs`) keeps the lists of
    // at most `x_flow_window` chunks waiting for the reduction; a subject's translation still
    // lives only while its worker searches it.
    #[cfg(feature = "parallel")]
    let mut parallel_batch = std::collections::VecDeque::new();
    #[cfg(feature = "parallel")]
    let diagnose = crate::utils::threading::diagnostics_enabled();
    #[cfg(feature = "parallel")]
    let active = AtomicUsize::new(0);
    #[cfg(feature = "parallel")]
    let peak = AtomicUsize::new(0);
    #[cfg(feature = "parallel")]
    let mut workers = std::collections::BTreeSet::new();
    #[cfg(feature = "parallel")]
    let search = |x_oid: usize| {
        if diagnose {
            let now = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(now, Ordering::SeqCst);
        }
        let result = subject_core_oid(x_oid);
        if diagnose {
            active.fetch_sub(1, Ordering::SeqCst);
        }
        (result, rayon::current_thread_index())
    };
    // Graduated 2026-10-22 (was `LOSAT_X_TBNOIDNEXT`): every chunk handed to the threads when
    // they are free (`x_oid_flow.rs`); the lists come back below chunk by chunk in OID order.
    #[cfg(feature = "parallel")]
    crate::utils::xreach::x_reach!("LOSAT_X_TBNOIDNEXT", x_chunk_bases.is_some());
    // LOSAT_X_TBNFLOWHELP: the threads that wait in the flow run the pool's jobs; no flow, no wait.
    #[cfg(feature = "parallel")]
    if super::x_oid_flow::x_flow_help().is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_TBNFLOWHELP", x_chunk_bases.is_some());
    }
    // LOSAT_X_TBNLOOPFREE: the threads of the flow free the buffers they kept between subjects when
    // their loops end; no flow, nothing freed here.
    #[cfg(feature = "parallel")]
    if super::x_oid_flow::x_loop_free_on() {
        crate::utils::xreach::x_reach!("LOSAT_X_TBNLOOPFREE", x_chunk_bases.is_some());
    }
    #[cfg(feature = "parallel")]
    let x_flow = x_chunk_bases.map(|bases| {
        super::x_oid_flow::XOidFlow::new(
            x_all_chunks(subjects, bases),
            x_flow_window(pool.map_or(1, SearchPool::threads)),
            pool.map_or(1, SearchPool::threads).saturating_sub(1),
        )
    });
    #[cfg(feature = "parallel")]
    let x_flow_search = |x_oid: usize| search(x_oid);
    // The OID-ordered reduction below, as one closure so that it can run inside the scope whose
    // other threads search the chunks of the flow.
    let mut reduce_subjects = || -> Result<()> {
        for oid in 0..subjects.count() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905:
            #[cfg(feature = "parallel")]
            let (preliminary, linked) = if parallel_selected {
                if parallel_batch.is_empty() {
                    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:338-342
                    // Dispatch point: the next chunk in OID order from the flow (`x_oid_flow.rs`).
                    // The subjects are searched by the same call as before and their lists are
                    // taken in OID order. NCBI's thread searches the subjects of its chunk one
                    // after another; here they are a parallel range inside the chunk's task, so
                    // a thread that runs out of chunks takes part of a chunk another thread is
                    // still searching instead of waiting for all of it (a subject longer than a
                    // chunk's budget is a chunk of its own). The ordered collect keeps the
                    // subject order. Scheduling only.
                    let flow = x_flow.as_ref().expect("the flow of a selected pool");
                    let batch: Vec<_> = flow
                        .next_in_order(&x_flow_search)?
                        .context("missing ordered TBLASTN subject chunk")?;
                    if diagnose {
                        workers.extend(batch.iter().filter_map(|(_, worker)| *worker));
                    }
                    parallel_batch.extend(batch);
                }
                parallel_batch
                    .pop_front()
                    .context("missing ordered TBLASTN subject result")?
                    .0?
            } else {
                subject_core_oid(oid)?
            };
            #[cfg(not(feature = "parallel"))]
            let (preliminary, linked) = subject_core_oid(oid)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905:
            // The comparison-only trace retains pre-link HSP bytes by OID.
            if let Some(ref mut trace) = observer {
                trace.preliminary.push((oid, preliminary.clone()));
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1555:
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:289-319:
            // Query-indexed lists are consumed in
            // reverse-sorted order by BlastHSPStreamRead.
            sink.add(oid, &linked)?;
        }
        Ok(())
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-180
    // NCBI runs its search threads side by side until the
    // subjects run out. With the flow the other threads of the pool take chunks for the whole
    // search while this thread reduces (and searches when the next chunk is not ready); the flow
    // is closed when the reduction ends, fails or panics, so the scope always ends.
    #[cfg(feature = "parallel")]
    match (x_flow.as_ref(), pool) {
        (Some(flow), Some(pool)) => pool.install(|| {
            let reduced = rayon::scope(|scope| {
                for _ in 1..pool.threads() {
                    scope.spawn(|_| flow.work(&x_flow_search));
                }
                let _close = super::x_oid_flow::XCloseOnDrop(flow);
                reduce_subjects()
            });
            // LOSAT_X_TBNLOOPFREE: the owner's loop has ended too.
            super::x_oid_flow::x_loop_free();
            reduced
        })?,
        _ => reduce_subjects()?,
    }
    #[cfg(not(feature = "parallel"))]
    reduce_subjects()?;
    #[cfg(feature = "parallel")]
    if let Some(line) = x_flow.as_ref().and_then(|flow| flow.stats_line()) {
        eprintln!("{line}");
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1659-1668:
    //    argument can be shared between multiple threads */
    #[cfg(feature = "parallel")]
    if parallel_selected && diagnose {
        // The search is one entry into the pool; the line adds its chunks and the most finished
        // chunks that waited for the reduction at once.
        let x_flow_fields = x_flow.as_ref().map_or(String::new(), |flow| {
            format!(
                " flow_chunks={} flow_peak_waiting={}",
                flow.chunks(),
                flow.peak_waiting()
            )
        });
        crate::utils::threading::write_diagnostic_line(format!(
            "[losat-thread-activity] program=tblastn stage=subject_core work_items={} worker_slots={:?} peak_active={} pool_entries=1{}\n",
            subjects.count(), workers, peak.load(Ordering::SeqCst), x_flow_fields
        ));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:728-729,804-813,882-887:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:508-527,1070-1101:
// frame_offsets[context+1] = offset after the translated frame and sentinel.
// The final nonempty frame in the six-frame loop supplies stat_length.
fn initial_translated_stat_length(nt_length: usize) -> usize {
    match nt_length {
        0..=2 => nt_length,
        n => ((n - 2) / 3).max(1),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:141-144,1811-1926:
// Blast_HSPNew() uses calloc, so num remains zero without BLAST_LinkHsps.
// hsp->evalue = BLAST_SpougeStoE(score, kbp, sbp->gbp,
fn direct_local_evalues(
    input: &[(usize, super::search_gapped::GappedHsp)],
    query_lengths: &[i32],
    subject_stat_length: i32,
    gapped_params: &[KarlinParams],
    lengths: &[super::stage_d_stats::LocalContextLength],
    gumbel: Option<&BlastGumbelBlk>,
) -> LinkedHspList {
    let hsps: Vec<_> = input
        .iter()
        .enumerate()
        .map(|(source_index, &(context, hsp))| LinkedHsp {
            context,
            hsp,
            num: 0,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1870-1895:
            evalue: if let Some(gumbel) = gumbel {
                blast_spouge_stoe(
                    hsp.score,
                    &gapped_params[context],
                    gumbel,
                    query_lengths[context],
                    subject_stat_length,
                )
            } else {
                evalue_from_raw_score(
                    hsp.score,
                    &gapped_params[context],
                    lengths[context].eff_searchsp as f64,
                )
            },
            source_index,
        })
        .collect();
    let best_evalue = hsps
        .iter()
        .map(|h| h.evalue)
        .reduce(f64::min)
        .unwrap_or(0.0);
    LinkedHspList { hsps, best_evalue }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:200-251,675-721,
// 1644-1707,1763-1776; core/blast_hits.c:1811-1935,3243-3297:
fn finish_local_mode0_no_sum_stats(
    queries: &[Vec<u8>],
    subjects: SubjectLetters<'_>,
    profile: LocalStageDProfile<'_>,
    preliminary_lists: Vec<(usize, usize, LinkedHspList)>,
    query_lengths: &[i32],
    gapped_params: &[KarlinParams],
    gumbel: Option<&BlastGumbelBlk>,
    code: &GeneticCode,
    lengths: &[super::stage_d_stats::LocalContextLength],
    link: Option<&super::stage_d_stats::LocalLinkParameters>,
    scoring: LocalStageDScoring,
    mut observer: Option<&mut StageDBoundaryTrace>,
    pool: Option<&SearchPool<'_>>,
) -> Result<Vec<KappaResultHitList>> {
    let mut results: Vec<_> = (0..queries.len())
        .map(|_| KappaResultHitList::new(profile.max_target_seqs))
        .collect::<Result<_>>()?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-605;
    // c++/src/algo/blast/core/blast_encoding.c:110-121:
    // the traceback identity pass reads the unmasked amino-acid query;
    // lower-case FASTA residues represent the same protein letters.
    let query_nomask: Vec<_> = queries
        .iter()
        .map(|query| {
            query
                .iter()
                .copied()
                .map(|aa| aa_char_to_ncbistdaa(aa.to_ascii_uppercase()))
                .collect::<Vec<_>>()
        })
        .collect();
    // With a `-db` stream, the subjects of the lists, each read
    // when first aligned, held until their last list is done (LOSAT_X_TBNSTAGEDROP; with it off
    // until the traceback ends).
    let x_stage = super::x_db_stream::XStageSubjects::new(
        subjects,
        preliminary_lists.iter().map(|(_, oid, _)| *oid),
    );
    // The traceback of one list and what its result keeps: the posttraceback subject length, and
    // the HSPs with their statistics and report counts unless none is left. A function of the list
    // alone (its query, its subject and the search's parameters).
    let item = |context: usize,
                oid: usize,
                preliminary: &LinkedHspList|
     -> Result<(i32, Option<KappaResultList>)> {
        // The subject as the readers of windows take it: the record's letters, or its store (every
        // `-db` sequence; with EXPERIMENT LOSAT_X_TBNSTORE every `-subject` record), whose windows
        // the translations make when they read them (`TargetTranslation`), also for a list retried
        // after a fence, whose frames NCBI translates whole (windows of the whole translation,
        // `TargetTranslation::x_whole_window`). A `-db` sequence is
        // read once for the traceback (`x_stage`); EXPERIMENT (LOSAT_X_TBNSTAGEDROP): the stage
        // lets it go after its last list (`_done`, declared before `held`).
        let _done = x_stage.x_match(oid);
        let held = x_stage.held(oid)?;
        let subject = held.nt();
        let input: Vec<_> = preliminary.hsps.iter().map(|linked| linked.hsp).collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-536:
        // traceback uses the selected scoring matrix and gap costs.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-463
        // The final X-drop is at least the preliminary one (`-xdrop_gap` above
        // `-xdrop_gap_final`).
        let final_xdrop = local_extension_final_xdrop(
            scoring.gap_xdrop_bits,
            scoring.final_xdrop_bits,
            gapped_params[context].lambda,
        )?;
        // Dispatch point: the scores of the batch's PSSM query (`-in_pssm`, the same for every
        // list), else the matrix's.
        let (traced, stat_length) = match profile.pssm {
            None => full_translation_traceback_with_matrix_and_events_with_mask_mode_owned(
                &queries[context],
                subject,
                profile.genetic_code,
                &input,
                scoring.matrix,
                scoring.gap_open,
                scoring.gap_extend,
                final_xdrop,
                0.0,
                profile.seg,
                profile.soft_masking,
                profile.mask_lowercase,
                None,
                None,
                None,
                None,
            )?,
            Some(psi) => super::search_gapped::x_pssm_traceback_owned(
                &psi.pssm,
                &queries[context],
                subject,
                profile.genetic_code,
                &input,
                scoring.matrix,
                scoring.gap_open,
                scoring.gap_extend,
                final_xdrop,
                profile.seg,
                profile.soft_masking,
                profile.mask_lowercase,
            )?,
        };
        if traced.is_empty() {
            return Ok((stat_length, None));
        }
        let input: Vec<_> = traced.iter().map(|row| (context, row.hsp)).collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:213-250,717-719:
        // else Blast_HSPListGetEvalues(..., stat_length, ...);
        let mut post = if let Some(link) = link {
            let mut post_link = *link;
            post_link.cutoff_small_gap = 0;
            link_preliminary_hsps(
                &input,
                query_lengths,
                lengths,
                stat_length,
                gapped_params,
                gumbel,
                &post_link,
            )?
        } else {
            direct_local_evalues(
                &input,
                query_lengths,
                stat_length,
                gapped_params,
                lengths,
                gumbel,
            )
        };
        reap_by_evalue(&mut post, profile.expect_value);
        if post.hsps.is_empty() {
            return Ok((stat_length, None));
        }
        let mut target = TargetTranslation::of(subject, code);
        let mut owned: Vec<Option<TracebackOwnedHsp>> = traced.into_iter().map(Some).collect();
        let payloads = post
            .hsps
            .iter()
            .map(|hsp| {
                let original = owned[hsp.source_index]
                    .take()
                    .context("TBLASTN ordinary traceback source owner reused")?;
                let converted = ConvertedKappaHsp {
                    context,
                    hsp: original.hsp,
                    edit_script: original.edit_script,
                    matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
                };
                let (ident, positives, length, mismatches, gap_opens, gap_letters) =
                    postredo_converted_stats_with_matrix(
                        &converted,
                        &query_nomask[context],
                        &mut target,
                        scoring.matrix,
                    )?;
                Ok(KappaHspPayload {
                    bit_score: (hsp.hsp.score as f64 * gapped_params[context].lambda
                        - gapped_params[context].k.glibc_ln())
                        / std::f64::consts::LN_2,
                    num_ident: i32::try_from(ident)?,
                    num_positives: positives,
                    report_num_ident: ident,
                    report_num_positives: positives,
                    report_mismatches: i64::try_from(mismatches)?,
                    align_length: length,
                    mismatches,
                    gap_opens,
                    gap_letters,
                    edit_script: converted.edit_script,
                    matrix_adjust_rule: converted.matrix_adjust_rule,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok((
            stat_length,
            Some(KappaResultList {
                oid: i32::try_from(oid)?,
                hsps: post,
                payloads,
            }),
        ))
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1518-1521
    // NCBI traces the batches of its HSP stream (the lists of one subject) on its threads
    // (`schedule(guided)`), each into the results of its thread, which are joined afterward
    // (SThreadLocalDataArrayConsolidateResults). Here the lists of a batch are traced on the search
    // pool and their results put in the hit lists in the order of the lists, as the serial loop
    // puts them; a list's result is a function of the list alone.
    #[cfg(feature = "parallel")]
    let parallel = pool.is_some_and(SearchPool::enabled) && preliminary_lists.len() > 1;
    #[cfg(not(feature = "parallel"))]
    let parallel = {
        let _ = pool;
        false
    };
    // No NCBI counterpart: stage report for the thread-use statistics; it does not change any value
    // NCBI computes.
    crate::utils::threading::report_stage(
        "tblastn",
        "traceback",
        preliminary_lists.len(),
        parallel,
    );
    let mut keep = |context: usize, oid: usize, (stat_length, result)| -> Result<()> {
        if let Some(ref mut trace) = observer {
            trace
                .posttraceback_lengths
                .push((context, oid, stat_length));
        }
        if let Some(result) = result {
            results[context].update(result)?;
        }
        Ok(())
    };
    if parallel {
        #[cfg(feature = "parallel")]
        {
            let pool = pool.expect("parallel pool selected");
            let batch_size = pool.threads().saturating_mul(32).max(64);
            for batch in preliminary_lists.chunks(batch_size) {
                let traced: Vec<_> = pool.install(|| {
                    batch
                        .par_iter()
                        .map(|(context, oid, list)| item(*context, *oid, list))
                        .collect()
                });
                for ((context, oid, _), result) in batch.iter().zip(traced) {
                    keep(*context, *oid, result?)?;
                }
            }
        }
    } else {
        for (context, oid, list) in &preliminary_lists {
            keep(*context, *oid, item(*context, *oid, list)?)?;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1776:
    // Otherwise the lists keep the order of their insertion above, by OID, then the heap's
    // (`KappaResultHitList::update`): NCBI's one-thread traceback.
    if profile.subject_set.sorts_results() {
        for result in &mut results {
            result.sort_by_evalue_for_positive_totlen();
        }
    }
    Ok(results)
}

#[cfg(test)]
#[path = "stage_d_pipeline_tests.rs"]
mod tests;
