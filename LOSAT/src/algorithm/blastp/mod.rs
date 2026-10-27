//! BLASTP algorithm module
//!
//! This module implements BLASTP (protein vs protein search) in pure Rust,
//! using NCBI BLAST as the behavioral reference.

mod alignment;
pub(crate) mod archive;
pub mod args;
pub mod blast_engine;
// NCBI c++/src/algo/blast/core/blast_setup.c:382-385:
// if (Blast_QueryIsProtein(program_number)) BLAST_ScoreSetAmbigRes(sbp, 'X');
// BLASTP and TBLASTN both encode an untranslated protein query in NCBISTDAA.
pub(crate) mod encoding;
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:217-228
// status = s_BlastAaWordFinder_TwoHit(subject, query, ..., matrix, ...);
// The same protein two-hit extension primitive is used by BLASTP and TBLASTN.
pub(crate) mod extension;
// The ungapped X-drop walks of aa_ungapped.c:804-883 with a byte score table, same results.
pub(crate) mod x_ungapped_fast;
// EXPERIMENT (LOSAT_X_BPUNGAPSENT): blastp's two-hit extension for word size 3 inlined into the
// fused loop's candidate closure, walks stopped by the sentinels, same results.
pub(crate) mod x_ungapped_sent;
// NCBI c++/src/algo/blast/core/blast_gapalign.c:3936-3953:
// max_offset = BlastGetStartForGappedAlignment(...);
// status = s_BlastProtGappedAlignment(...);
// BLASTP and TBLASTN use the same protein gapped-score path.
mod blast_xml;
pub(crate) mod gapalign;
mod hsp;
mod query_split;
mod seqalign;
// EXPERIMENT (LOSAT_X_BPQBLOCK): the queries searched in blocks of whole query batches.
mod x_query_blocks;
// Redone HSP lists without spare room.
mod x_redo_fit;
// EXPERIMENT (LOSAT_X_BPREDOPAR): the composition-based redo with jobs
// of up to four matches on the pool, replayed in stream order. Built where the parallel redo paths
// are.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod x_redo_par;
// The redo of a `-db` search reads each match's subject from the database.
mod x_redo_subjects;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:102-3691
// ```c
// static void s_HSPListNormalizeScores(...);
// static void s_HitlistReapContained(...);
// static int s_HitlistEvaluateAndPurge(...);
// static void s_ComputeNumIdentities(...);
// ```
mod kappa;
// The subject scan and the two-hit test in one loop.
pub(crate) mod x_seed_fused;
// The same loop for the compressed lookup table of blastp-fast.
mod x_seed_fused_compressed;
// With LOSAT_X_BPCLUT: the same loop on the flat compressed table.
mod x_seed_fused_clut;
// Preliminary hit lists without spare room, each freed after its query's redo.
mod x_prelim_free;
// Preliminary hit lists with a 64-byte element until the redo.
mod x_prelim_compact;

/// EXPERIMENT (LOSAT_X_BPHITBF): the `[X_STATS] BPHITBF` line at exit.
// No NCBI counterpart: it prints a counter and changes no value NCBI computes.
pub fn x_bp_hit_bf_print_summary() {
    x_seed_fused::x_bp_hit_bf_print_summary();
}

/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the `[X_STATS] BPLKPAHEAD` line at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub fn x_bp_lkp_ahead_print_summary() {
    blast_engine::x_lookup_ahead::x_bp_lkp_ahead_print_summary();
}

pub use args::BlastpArgs;
pub use blast_engine::run;
// EXPERIMENT (LOSAT_X_BPGAP): the summary line `main` prints at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub use gapalign::x_restricted_gap_print_summary;
