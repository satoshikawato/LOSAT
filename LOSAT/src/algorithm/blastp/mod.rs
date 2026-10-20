//! BLASTP algorithm module
//!
//! This module implements BLASTP (protein vs protein search) in pure Rust,
//! using NCBI BLAST as the behavioral reference.

mod alignment;
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
// EXPERIMENT (LOSAT_X_BPUNGAP / LOSAT_X_BPUNGAPSHADOW): the ungapped X-drop
// walks of aa_ungapped.c:804-883 with a byte score table, same results (`pub` for the summary line
// of `main`).
pub mod x_ungapped_fast;
// NCBI c++/src/algo/blast/core/blast_gapalign.c:3936-3953:
// max_offset = BlastGetStartForGappedAlignment(...);
// status = s_BlastProtGappedAlignment(...);
// BLASTP and TBLASTN use the same protein gapped-score path.
pub(crate) mod gapalign;
mod hsp;
mod query_split;
mod seqalign;
// EXPERIMENT (LOSAT_X_BPQBLOCK): the queries searched in blocks of whole query batches.
mod x_query_blocks;
// EXPERIMENT (LOSAT_X_BPREDOFIT): redone HSP lists without spare room.
mod x_redo_fit;
// EXPERIMENT (LOSAT_X_BPREDOPAR / LOSAT_X_BPREDOPARSHADOW): the composition-based redo with jobs
// of up to four matches on the pool, replayed in stream order. Built where the parallel redo paths
// are.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod x_redo_par;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:102-3691
// ```c
// static void s_HSPListNormalizeScores(...);
// static void s_HitlistReapContained(...);
// static int s_HitlistEvaluateAndPurge(...);
// static void s_ComputeNumIdentities(...);
// ```
mod kappa;
// EXPERIMENT (LOSAT_X_BPSEED): the subject scan and the two-hit test in one loop.
pub(crate) mod x_seed_fused;
// EXPERIMENT (LOSAT_X_BPSEED): the same loop for the compressed lookup table of blastp-fast.
mod x_seed_fused_compressed;
// EXPERIMENT (LOSAT_X_BPSEED with LOSAT_X_BPCLUT): the same loop on the flat compressed table.
mod x_seed_fused_clut;
// EXPERIMENT (LOSAT_X_BPPFREE): preliminary hit lists without spare room, each freed after its
// query's redo.
mod x_prelim_free;
// EXPERIMENT (LOSAT_X_BPPCOMPACT): preliminary hit lists with a 64-byte element until the redo.
mod x_prelim_compact;

/// EXPERIMENT (LOSAT_X_BPSEEDSHADOW): the shadow mode's line at exit.
pub fn x_bp_seed_print_summary() {
    x_seed_fused::x_bp_seed_print_summary();
}

/// EXPERIMENT (LOSAT_X_BPLKPAHEADSHADOW): the shadow mode's line at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub fn x_bp_lkp_ahead_print_summary() {
    blast_engine::x_lookup_ahead::x_bp_lkp_ahead_print_summary();
}

/// EXPERIMENT (LOSAT_X_BPREDOPARSHADOW): the shadow mode's line at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub fn x_bp_redo_par_print_summary() {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    x_redo_par::x_bp_redo_par_print_summary();
}

/// EXPERIMENT (LOSAT_X_BPSUBJRUNSHADOW): the shadow mode's line at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub fn x_bp_subj_run_print_summary() {
    blast_engine::x_bp_subj_run_print_summary();
}

/// EXPERIMENT (LOSAT_X_BPDBCUTSHADOW): the shadow mode's line at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub fn x_bp_db_cut_print_summary() {
    blast_engine::x_bp_db_cut_print_summary();
}

/// EXPERIMENT (LOSAT_X_BPSUBJMETASHADOW): the shadow mode's line at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub fn x_bp_subj_meta_print_summary() {
    blast_engine::x_bp_subj_meta_print_summary();
}

/// EXPERIMENT (LOSAT_X_BPREDOPARSHADOW): a search's result against the new path's when the
/// one-thread redo of the search ended with an error (it returns before the shadow compares the
/// hit lists).
// No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
fn x_bp_redo_par_shadow_error(result: &anyhow::Result<()>) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    x_redo_par::x_redo_par_shadow_error(result);
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = result;
}

pub use args::BlastpArgs;
pub use blast_engine::run;
// EXPERIMENT (LOSAT_X_BPGAP / LOSAT_X_BPGAPSHADOW): the summary lines `main` prints at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub use gapalign::x_restricted_gap_print_summary;
// EXPERIMENT (LOSAT_X_BPPCOMPACT=shadow): the summary line `main` prints at exit.
pub use x_prelim_compact::x_bp_prelim_compact_print_summary;
