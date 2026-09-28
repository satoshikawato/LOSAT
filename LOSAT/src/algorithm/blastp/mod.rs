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
// NCBI c++/src/algo/blast/core/blast_gapalign.c:3936-3953:
// max_offset = BlastGetStartForGappedAlignment(...);
// status = s_BlastProtGappedAlignment(...);
// BLASTP and TBLASTN use the same protein gapped-score path.
pub(crate) mod gapalign;
mod hsp;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:102-3691
// ```c
// static void s_HSPListNormalizeScores(...);
// static void s_HitlistReapContained(...);
// static int s_HitlistEvaluateAndPurge(...);
// static void s_ComputeNumIdentities(...);
// ```
mod kappa;

pub use args::BlastpArgs;
pub use blast_engine::run;
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
// ```c
// SetupQueries_OMF(IBlastQuerySource& queries,
//                  BlastQueryInfo* qinfo,
//                  BLAST_SequenceBlk** seqblk,
//                  EBlastProgramType prog,
//                  ...)
// ```
#[cfg(target_arch = "wasm32")]
pub use blast_engine::run_web_pair;
#[cfg(target_arch = "wasm32")]
pub use blast_engine::run_web_pair_records;
