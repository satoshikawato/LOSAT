pub mod dust;
pub mod genetic_code;
pub mod matrix;
pub mod seg;
mod seg_lnfact;
// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
// (*thread)->Run(); (*thread)->Join(&result);
pub mod threading;
pub mod winmask;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1533-1537
// tid = omp_get_thread_num(); ... gap_align = thread_data->tld[tid]->gap_align;
pub(crate) mod scratch_pool;
pub(crate) mod x87;
pub(crate) mod x_frame;
pub mod x_logclone;
pub(crate) mod x_mapped;
pub mod x_seg_trim_memo;
pub(crate) mod xahead;
pub mod xdiag;
pub(crate) mod xdrop_simd;
pub(crate) mod xenv;
pub mod xreach;
pub mod xstats;
pub(crate) mod xswitch;
pub mod xtiming;
