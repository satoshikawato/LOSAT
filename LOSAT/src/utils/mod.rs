pub mod dust;
pub mod genetic_code;
pub mod matrix;
pub mod seg;
mod seg_lnfact;
// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
// (*thread)->Run(); (*thread)->Join(&result);
pub mod threading;
// The global allocator of the threaded Wasm builds; its tests run on every target.
#[cfg(any(test, all(target_arch = "wasm32", feature = "wasm-threads")))]
mod wasm_alloc;
pub mod winmask;
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
